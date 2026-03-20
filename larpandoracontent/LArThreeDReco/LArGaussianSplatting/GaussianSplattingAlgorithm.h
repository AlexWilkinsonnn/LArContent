/**
 *  @file   larpandoracontent/LArThreeDReco/LArGaussianSplatting/GaussianSplattingAlgorithm.h
 *
 *  @brief  Header file for the Gaussian splatting 2D->3D matching algorithm.
 *
 *  Fits 3D Gaussian primitives to 2D clusters across U, V, W wire-plane views
 *  using gradient descent, then computes asymmetric cluster-cluster coverage
 *  scores via shared Gaussian responsibilities.
 *
 *  $Log: $
 */
#ifndef LAR_GAUSSIAN_SPLATTING_ALGORITHM_H
#define LAR_GAUSSIAN_SPLATTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

#include <Eigen/Dense>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace lar_content
{

/**
 *  @brief  GaussianSplattingAlgorithm class
 *
 *  For each connected X-overlap region this algorithm:
 *    1. Initialises one 3D Gaussian per 2D cluster
 *    2. Fits all Gaussians simultaneously to hit data across all three views
 *       using Adam gradient descent
 *    3. Prunes low-opacity Gaussians
 *    4. Computes asymmetric coverage scores F(c_i->c_j) and F(c_j->c_i)
 *       for every cross-view cluster pair, which encode how much of each
 *       cluster's Gaussian content is shared with the other
 */
class GaussianSplattingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    GaussianSplattingAlgorithm();

    /**
     *  @brief  Asymmetric coverage scores between two clusters
     */
    struct CoverageScores
    {
        float m_forwardCoverage;  ///< F(c_i -> c_j): fraction of c_i's Gaussian content in c_j
        float m_reverseCoverage;  ///< F(c_j -> c_i): fraction of c_j's Gaussian content in c_i
    };

    /**
     *  @brief  Key for the coverage score map: ordered pair of cluster pointers
     */
    typedef std::pair<const pandora::Cluster *, const pandora::Cluster *> ClusterPair;
    typedef std::map<ClusterPair, CoverageScores> CoverageScoreMap;

    /**
     *  @brief  Get coverage scores computed during the last Run() call
     */
    const CoverageScoreMap &GetCoverageScores() const;

private:
    // -------------------------------------------------------------------------
    //  Internal data structures
    // -------------------------------------------------------------------------

    /**
     *  @brief  2D hit representation for the optimisation
     */
    struct Hit2D
    {
        Eigen::Vector2d m_position; ///< (x, wire_coord)
        double m_sigmaX;            ///< Hit width in drift direction (1-sigma)
        double m_sigmaWire;         ///< Hit width in wire direction (1-sigma)
    };

    /**
     *  @brief  All data belonging to one 2D cluster in one view
     */
    struct ClusterData
    {
        const pandora::Cluster *m_pCluster; ///< Pointer to the Pandora cluster
        pandora::HitType m_hitType;         ///< Which wire-plane view
        std::vector<Hit2D> m_hits;          ///< Extracted hit data
        float m_xMin;                       ///< Cluster X extent minimum
        float m_xMax;                       ///< Cluster X extent maximum
        float m_wireMin;                    ///< Cluster wire-coord extent minimum
        float m_wireMax;                    ///< Cluster wire-coord extent maximum
    };

    /**
     *  @brief  A group of clusters from all views sharing an X-overlap region
     */
    struct ClusterGroup
    {
        std::vector<ClusterData> m_clusters; ///< All clusters in the overlap region
        float m_xOverlapMin;                 ///< Common X-overlap minimum
        float m_xOverlapMax;                 ///< Common X-overlap maximum
    };

    /**
     *  @brief  A 3D Gaussian primitive with diagonal covariance
     *
     *  Covariance is parametrised as diag(exp(2*lx), exp(2*ly), exp(2*lz))
     *  so that the unconstrained log-sigma parameters can be freely optimised.
     */
    struct Gaussian3D
    {
        Eigen::Vector3d m_mean;     ///< μ = (x, y, z)
        Eigen::Vector3d m_logSigma; ///< log(σ_x, σ_y, σ_z); sigma = exp(logSigma)
        double m_opacity;           ///< α ∈ (0, 1]; pruned when below threshold

        /**
         *  @brief  Get diagonal covariance matrix Σ = diag(σ_x², σ_y², σ_z²)
         */
        Eigen::Matrix3d GetCovariance() const;

        /**
         *  @brief  Get Gaussian volume (2π)^{3/2} * |Σ|^{1/2} = (2π)^{3/2} * σ_x σ_y σ_z
         */
        double GetVolume() const;
    };

    /**
     *  @brief  2×3 wire-plane projection matrix P_k
     *
     *  Maps 3D position (x,y,z) -> 2D (x, wire_coord) for a given view.
     *  Computed numerically from LArGeometryHelper::ProjectPosition.
     */
    struct ProjectionMatrix
    {
        Eigen::Matrix<double, 2, 3> m_P; ///< The 2x3 matrix
    };

    /**
     *  @brief  Per-parameter Adam optimiser state for one Gaussian
     */
    struct AdamState
    {
        Eigen::VectorXd m_firstMoment;  ///< First moment estimate (mean)
        Eigen::VectorXd m_secondMoment; ///< Second moment estimate (uncentred variance)
        int m_step;                     ///< Number of updates applied
    };

    // -------------------------------------------------------------------------
    //  Algorithm interface
    // -------------------------------------------------------------------------
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // -------------------------------------------------------------------------
    //  Pipeline methods
    // -------------------------------------------------------------------------

    /**
     *  @brief  Retrieve and filter clusters from the configured list names
     */
    void GetClusterData(const std::string &listName, pandora::HitType hitType, std::vector<ClusterData> &clusterData) const;

    /**
     *  @brief  Extract 2D hit data from a cluster
     */
    void ExtractHits(const pandora::Cluster *const pCluster, pandora::HitType hitType, ClusterData &clusterData) const;

    /**
     *  @brief  Group clusters from all three views into connected X-overlap regions
     */
    void GroupClustersByOverlap(const std::vector<ClusterData> &allClusters, std::vector<ClusterGroup> &groups) const;

    /**
     *  @brief  Compute the 2x3 projection matrix for a given wire-plane view
     */
    ProjectionMatrix ComputeProjectionMatrix(pandora::HitType hitType) const;

    /**
     *  @brief  Initialise one Gaussian per cluster, placed at the cluster centroid
     *          with the ambiguous coordinate set to the midpoint of the overlap region
     */
    void InitialiseGaussians(const ClusterGroup &group,
        const std::map<pandora::HitType, ProjectionMatrix> &projections, std::vector<Gaussian3D> &gaussians) const;

    /**
     *  @brief  Run Adam gradient descent to fit all Gaussians to the cluster hit data
     */
    void RunGradientDescent(const ClusterGroup &group,
        const std::map<pandora::HitType, ProjectionMatrix> &projections, std::vector<Gaussian3D> &gaussians) const;

    /**
     *  @brief  Compute the total loss and per-Gaussian gradients
     *
     *  Loss = -Σ_views Σ_hits Σ_gaussians α_g * N(hit; P_k μ_g, P_k Σ_g P_k^T + Σ_hit)
     *
     *  @param  group        cluster group containing all hits
     *  @param  projections  per-view projection matrices
     *  @param  gaussians    current Gaussian parameters
     *  @param  gradMean     output gradient w.r.t. mean for each Gaussian
     *  @param  gradLogSigma output gradient w.r.t. log-sigma for each Gaussian
     *  @param  gradOpacity  output gradient w.r.t. opacity for each Gaussian
     *
     *  @return total loss value
     */
    double ComputeLossAndGradients(const ClusterGroup &group, const std::map<pandora::HitType, ProjectionMatrix> &projections,
        const std::vector<Gaussian3D> &gaussians, std::vector<Eigen::Vector3d> &gradMean,
        std::vector<Eigen::Vector3d> &gradLogSigma, std::vector<double> &gradOpacity) const;

    /**
     *  @brief  Apply one Adam update step to all Gaussian parameters
     *
     *  @param  gaussians    Gaussian parameters (modified in place)
     *  @param  adamStates   per-Gaussian Adam moment estimates (modified in place)
     *  @param  gradMean     gradient w.r.t. mean
     *  @param  gradLogSigma gradient w.r.t. log-sigma
     *  @param  gradOpacity  gradient w.r.t. opacity
     */
    void ApplyAdamStep(std::vector<Gaussian3D> &gaussians, std::vector<AdamState> &adamStates,
        const std::vector<Eigen::Vector3d> &gradMean, const std::vector<Eigen::Vector3d> &gradLogSigma,
        const std::vector<double> &gradOpacity) const;

    /**
     *  @brief  Remove Gaussians whose opacity has fallen below m_opacityPruneThreshold
     */
    void PruneGaussians(std::vector<Gaussian3D> &gaussians) const;

    /**
     *  @brief  Visualise the 2D projections of fitted Gaussians for all three views.
     *
     *  For each Gaussian and each view, projects the 3D mean and covariance into 2D,
     *  then draws two lines along the major and minor axes of the projected ellipse
     *  (each line extends ±1σ from the projected centre).
     *  Only compiled when MONITORING is defined.
     *
     *  @param  gaussians    fitted (and pruned) Gaussian primitives
     *  @param  projections  per-view projection matrices
     */
    void VisualiseGaussians(const std::vector<Gaussian3D> &gaussians,
        const std::map<pandora::HitType, ProjectionMatrix> &projections) const;

    /**
     *  @brief  Compute the responsibility r(g, c_k): the fraction of Gaussian g's
     *          projected probability mass (in view k) that falls within cluster c_k's
     *          bounding box
     *
     *  Uses the product of 1D Gaussian CDFs over the cluster's [xMin,xMax] x [wireMin,wireMax].
     *
     *  @param  g          the 3D Gaussian
     *  @param  cd         the cluster data (defines bounding box and view)
     *  @param  projection the projection matrix for the cluster's view
     *
     *  @return responsibility r(g, c_k) in [0,1]
     */
    double ComputeResponsibility(
        const Gaussian3D &g, const ClusterData &cd, const ProjectionMatrix &projection) const;

    /**
     *  @brief  Compute asymmetric coverage scores for all cross-view cluster pairs
     *          in the group and insert into scoreMap
     *
     *  F(c_i -> c_j) = S(c_i, c_j) / S_self(c_i)  where
     *  S(c_i, c_j)  = Σ_g r(g,c_i) * r(g,c_j) * α_g * V_g   (volume-weighted shared affinity)
     *  S_self(c_i)  = Σ_g r(g,c_i)^2 * α_g * V_g
     *
     *  @param  group       the cluster group
     *  @param  projections per-view projection matrices
     *  @param  gaussians   fitted Gaussian primitives
     *  @param  scoreMap    output map (modified in place)
     */
    void ComputeCoverageScores(const ClusterGroup &group, const std::map<pandora::HitType, ProjectionMatrix> &projections,
        const std::vector<Gaussian3D> &gaussians, CoverageScoreMap &scoreMap) const;

    // -------------------------------------------------------------------------
    //  Configuration parameters
    // -------------------------------------------------------------------------
    std::string m_clusterListNameU; ///< U-view cluster list name
    std::string m_clusterListNameV; ///< V-view cluster list name
    std::string m_clusterListNameW; ///< W-view cluster list name

    unsigned int m_minClusterHits;   ///< Minimum number of hits in a cluster to consider
    float m_minClusterLength;        ///< Minimum cluster length to consider (cm)
    float m_xOverlapThreshold;       ///< X-gap tolerance: intervals closer than this are merged into one group (cm)
    float m_matchThreshold;          ///< Coverage score threshold for cluster matching decisions

    unsigned int m_nOptimisationSteps;  ///< Number of Adam gradient descent steps
    float m_learningRate;               ///< Adam learning rate
    float m_adamBeta1;                  ///< Adam first-moment decay rate
    float m_adamBeta2;                  ///< Adam second-moment decay rate
    float m_adamEpsilon;                ///< Adam numerical stability constant

    float m_initialSigmaXY;            ///< Initial Gaussian σ in x and y (cm)
    float m_initialSigmaWire;          ///< Initial Gaussian σ in wire direction (cm)
    float m_opacityPruneThreshold;     ///< Prune Gaussians with opacity below this

    // -------------------------------------------------------------------------
    //  Outputs
    // -------------------------------------------------------------------------
    mutable CoverageScoreMap m_coverageScores; ///< Coverage scores from latest Run()
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const GaussianSplattingAlgorithm::CoverageScoreMap &GaussianSplattingAlgorithm::GetCoverageScores() const
{
    return m_coverageScores;
}

} // namespace lar_content

#endif // #ifndef LAR_GAUSSIAN_SPLATTING_ALGORITHM_H
