/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the deep learning shower growing algorithm
 *
 *  $Log: $
 */
#ifndef LAR_DL_SHOWER_GROWING_ALGORITHM
#define LAR_DL_SHOWER_GROWING_ALGORITHM 1

#include <torch/torch.h>

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

// using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DLShowerGrowingAlgorithm class
 */
class DLShowerGrowingAlgorithm : public pandora::Algorithm
{
private:
    struct HitFeatures
    {
        HitFeatures();

        float m_xRel;
        float m_zRel;
        float m_rRel;
        float m_cosThetaRel;
        float m_sinThetaRel;
        float m_distToXGap;
        float m_xWidth; 
        float m_energy;
    };

    struct ClusterGroup
    {
        ClusterGroup();

        void insert(const pandora::Cluster *pCluster);
        const pandora::Cluster *GetRepresentativeCluster() const { return m_representativeCluster; }
        const std::unordered_set<const pandora::Cluster *> &GetClusters() const { return m_clusters; }
        size_t size() const { return m_clusters.size(); }
        bool empty() const { return m_clusters.empty(); }
        auto begin() const { return m_clusters.begin(); }
        auto end() const { return m_clusters.end(); }

        std::unordered_set<const pandora::Cluster *> m_clusters;
        const pandora::Cluster *m_representativeCluster;
    };

public:
    /**
     *  @brief Default constructor
     */
    DLShowerGrowingAlgorithm();

    /**
     *  @brief Default destructor
     */
    virtual ~DLShowerGrowingAlgorithm();

private:
    typedef std::map<const pandora::Cluster *const, std::map<const pandora::Cluster *const, float>> SimilarityMatrix;
    typedef std::map<const pandora::Cluster *const, std::vector<const pandora::Cluster *>> AdjacencyLists;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create the input files for training the network
     */ 
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief Do an inference
     */ 
    pandora::StatusCode Infer();

    /* Start training sample preparation methods */

    const pandora::MCParticle* GetMainMC(
        const pandora::CaloHit *const pCaloHit,
        std::map<const pandora::MCParticle *const, const pandora::MCParticle *const> &mcFoldTo) const;

    const pandora::MCParticle* FoldMCTo(const pandora::MCParticle *const pMC) const;

    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendentElectrons) const;

    const pandora::MCParticle* FoldPotentialDeltaRayTo(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pMC) const;

    bool IsEM(const pandora::MCParticle *const pMC) const;

    /* End training sample preparation methods */

    /* Start general helpers */

    std::map<pandora::HitType, pandora::CartesianVector> Get2DVertices() const;

    pandora::ClusterList GetAllClusters() const;

    pandora::StatusCode GetClusters(const std::string clusterListName, pandora::ClusterList &clusterList) const;

    pandora::StatusCode CalculateHitFeatures(
        const pandora::CaloHit *const pCaloHit, const pandora::CartesianVector vtxPos, HitFeatures &hitFeatures) const;

    /* End general helpers */

    /* Start inference methods */

    pandora::StatusCode PredictClusterSimilarityMatrix(
        const pandora::ClusterList &clusterList,
        const pandora::HitType view,
        const pandora::CartesianVector &vtxPos,
        SimilarityMatrix &clusterSimMat);

    pandora::StatusCode MakeClusterTensor(
        const std::vector<HitFeatures> &clusterFeatures, const pandora::HitType view, torch::Tensor &tensorCluster) const;

    pandora::StatusCode PopulateClusterSimilarityMatrix(
        const torch::Tensor &tensorSimMat, const pandora::ClusterList &clusterList, SimilarityMatrix &clusterSimMat) const;

    pandora::StatusCode ClusterFromSimilarity(
        const SimilarityMatrix &clusterSimMat, const float similarityThreshold, std::vector<ClusterGroup> &clusterGroups) const;

    pandora::StatusCode PopulateAdjacencyLists(
        const SimilarityMatrix &simMat,
        const float similarityThreshold,
        AdjacencyLists &coreClusterAdjLists,
        AdjacencyLists &accClusterAdjLists) const;

    pandora::StatusCode CalculateConnectedGroups(const AdjacencyLists &clusterAdjLists, std::vector<ClusterGroup> &clusterGroups) const;

    pandora::StatusCode AddClustersToGroups(
        const SimilarityMatrix &clusterSimMat,
        AdjacencyLists &ungroupedClusterAdjLists,
        const float similarityThreshold,
        std::vector<ClusterGroup> &clusterGroups) const;

    pandora::StatusCode PopulateSuperClusterSimilarityMatrix(
        const std::vector<ClusterGroup> &clusterGroups, const SimilarityMatrix &clusterSimMat, SimilarityMatrix &superClusterSimMat) const;

    pandora::StatusCode ExpandSuperClusterGroups(
        const std::vector<ClusterGroup> &superClusterGroups,
        const std::vector<ClusterGroup> &clusterGroups,
        std::vector<ClusterGroup> &expandedSuperClusterGroups) const;

    pandora::StatusCode MergeGroups(const std::vector<ClusterGroup> &clusterGroups, const std::string &listNmae) const;

    /* End inference methods */

    /* Start shared mutable members */

    LArDLHelper::TorchModel m_modelEncoder; ///< TorchScript model for encoding hits in a cluster, set by "ModelEncoderFileName"
    LArDLHelper::TorchModel m_modelAttn;    ///< TorchScript model for attention over encoded clusters in a view, set by "ModelAttnFileName"
    LArDLHelper::TorchModel m_modelSim;     ///< TorchScript model for pairwise similarities over clusters after attention, set by "ModelSimFileName"
    float m_polarRScaleFactor;              ///< Scale factor for polar r coordinate input features
    float m_cartesianXScaleFactor;          ///< Scale factor for cartesian x coordinate input features
    float m_cartesianZScaleFactor;          ///< Scale factor for cartesian z coordinate input features
    std::set<double> m_detectorXGaps;       ///< X coordinates where gaps in X direction start/end

    /* End shared mutable members */

    /* Start hardcoded members */

    std::map<pandora::HitType, float> m_deltaRayLengthThresholdSquared; ///< Threshold for defining small delta rays that will be folded to the parent particle
    float m_deltaRayParentWeightThreshold;                              ///< Threshold for weight contribution of parent particle for it take the delta ray's hit
    int m_hitFeatureDim;                                                ///< Feature dimensions of each hit

    /* End hardcoded members */

    /* Start configurable via xml members */

    bool m_trainingMode;                       ///< Training mode
    std::string m_trainingTreeName;            ///< Tree name for training data output
    std::string m_trainingFileName;            ///< File name for training data output
    pandora::StringVector m_clusterListNames;  ///< Names of cluster lists
    std::string m_vertexListName;              ///< Name of vertex list
    float m_similarityThreshold;               ///< Threshold value on similarity for clusters to be connected
    float m_similarityThresholdBeta;           ///< Scaling factor for second clustering pass
    unsigned int m_accessoryClustersMaxHits;   ///< Clusters with this number of less hits are treated as accessory clusters during merging

    /* End configurable via xml members */
};

} // namespace lar_dl_content

#endif // LAR_DL_SHOWER_GROWING_ALGORITHM
