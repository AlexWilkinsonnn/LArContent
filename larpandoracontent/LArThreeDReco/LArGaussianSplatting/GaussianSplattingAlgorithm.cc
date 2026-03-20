/**
 *  @file   larpandoracontent/LArThreeDReco/LArGaussianSplatting/GaussianSplattingAlgorithm.cc
 *
 *  @brief  Implementation of the Gaussian splatting 2D->3D matching algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArGaussianSplatting/GaussianSplattingAlgorithm.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------
// Gaussian3D helpers
//------------------------------------------------------------------------------------------------------------------------------------------

Eigen::Matrix3d GaussianSplattingAlgorithm::Gaussian3D::GetCovariance() const
{
    const Eigen::Vector3d sigma = m_logSigma.array().exp();
    return (sigma.array() * sigma.array()).matrix().asDiagonal();
}

//------------------------------------------------------------------------------------------------------------------------------------------

double GaussianSplattingAlgorithm::Gaussian3D::GetVolume() const
{
    // (2π)^{3/2} * σ_x * σ_y * σ_z
    const Eigen::Vector3d sigma = m_logSigma.array().exp();
    return std::pow(2.0 * M_PI, 1.5) * sigma.prod();
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------------------------------------------------------------------

GaussianSplattingAlgorithm::GaussianSplattingAlgorithm() :
    m_minClusterHits(5),
    m_minClusterLength(0.f),
    m_xOverlapThreshold(0.5f),
    m_matchThreshold(0.5f),
    m_nOptimisationSteps(300),
    m_learningRate(0.05f),
    m_adamBeta1(0.9f),
    m_adamBeta2(0.999f),
    m_adamEpsilon(1e-8f),
    m_initialSigmaXY(1.f),
    m_initialSigmaWire(0.5f),
    m_opacityPruneThreshold(0.01f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Run
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GaussianSplattingAlgorithm::Run()
{
    m_coverageScores.clear();

    // Collect cluster data from all three views
    std::vector<ClusterData> allClusters;
    this->GetClusterData(m_clusterListNameU, TPC_VIEW_U, allClusters);
    this->GetClusterData(m_clusterListNameV, TPC_VIEW_V, allClusters);
    this->GetClusterData(m_clusterListNameW, TPC_VIEW_W, allClusters);

    if (allClusters.empty())
        return STATUS_CODE_SUCCESS;

    // Pre-compute projection matrices for each view
    std::map<HitType, ProjectionMatrix> projections;
    projections[TPC_VIEW_U] = this->ComputeProjectionMatrix(TPC_VIEW_U);
    projections[TPC_VIEW_V] = this->ComputeProjectionMatrix(TPC_VIEW_V);
    projections[TPC_VIEW_W] = this->ComputeProjectionMatrix(TPC_VIEW_W);

    // Group clusters by connected X-overlap regions
    std::vector<ClusterGroup> groups;
    this->GroupClustersByOverlap(allClusters, groups);
    std::cout << allClusters.size() << "\n";
    std::cout << groups.size() << "\n";
    std::set<const Cluster*> groupedClusterSet;
    ClusterList groupedClusterList;
    for (ClusterGroup group : groups)
    {
        std::cout << group.m_clusters.size() << ", ";
        for (ClusterData clusterData : group.m_clusters)
        {
            groupedClusterSet.insert(clusterData.m_pCluster);
            groupedClusterList.emplace_back(clusterData.m_pCluster);
        }
    }
    std::cout << "\n" << groupedClusterList.size() << ", " << groupedClusterSet.size() << "\n";

    // Process each group independently
    for (const ClusterGroup &group : groups)
    {
        std::map<HitType, ClusterList> viewToClusters;
        for (ClusterData clusterData : group.m_clusters)
            viewToClusters[clusterData.m_hitType].emplace_back(clusterData.m_pCluster);
        for (const auto &[view, clusterList] : viewToClusters)
        {
            const std::string name{"cluster group, view " + std::to_string(view)};
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, name.c_str(), AUTOITER));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        std::vector<Gaussian3D> gaussians;
        this->InitialiseGaussians(group, projections, gaussians);

        this->VisualiseGaussians(gaussians, projections);
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        if (gaussians.empty())
            continue;

        this->RunGradientDescent(group, projections, gaussians);
        this->PruneGaussians(gaussians);

        if (gaussians.empty())
            continue;

        this->VisualiseGaussians(gaussians, projections);
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        this->ComputeCoverageScores(group, projections, gaussians, m_coverageScores);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// GetClusterData
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::GetClusterData(
    const std::string &listName, const HitType hitType, std::vector<ClusterData> &clusterData) const
{
    const ClusterList *pClusterList{nullptr};

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pClusterList))
        return;

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (m_minClusterHits > 0 && pCluster->GetNCaloHits() < m_minClusterHits)
            continue;

        if (m_minClusterLength > 0.f && LArClusterHelper::GetLength(pCluster) < m_minClusterLength)
            continue;

        ClusterData cd;
        cd.m_pCluster = pCluster;
        cd.m_hitType = hitType;
        this->ExtractHits(pCluster, hitType, cd);

        if (cd.m_hits.empty())
            continue;

        clusterData.push_back(std::move(cd));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::ExtractHits(
    const Cluster *const pCluster, const HitType /*hitType*/, ClusterData &clusterData) const
{
    float xMin(std::numeric_limits<float>::max()), xMax(-std::numeric_limits<float>::max());
    float wireMin(std::numeric_limits<float>::max()), wireMax(-std::numeric_limits<float>::max());

    const OrderedCaloHitList &orderedList(pCluster->GetOrderedCaloHitList());

    for (const auto &layerEntry : orderedList)
    {
        for (const CaloHit *const pHit : *layerEntry.second)
        {
            const CartesianVector &pos(pHit->GetPositionVector());
            const float x(pos.GetX());
            const float wire(pos.GetZ()); // In 2D view: y=0, z = wire coordinate

            // CellSize0 = wire pitch; CellSize1 = drift sampling width
            const float sigmaX(pHit->GetCellSize1() * 0.5f);
            const float sigmaWire(pHit->GetCellSize0() * 0.5f);

            Hit2D hit;
            hit.m_position = Eigen::Vector2d(static_cast<double>(x), static_cast<double>(wire));
            hit.m_sigmaX = std::max(static_cast<double>(sigmaX), 0.01);
            hit.m_sigmaWire = std::max(static_cast<double>(sigmaWire), 0.01);
            clusterData.m_hits.push_back(hit);

            xMin = std::min(xMin, x);
            xMax = std::max(xMax, x);
            wireMin = std::min(wireMin, wire);
            wireMax = std::max(wireMax, wire);
        }
    }

    clusterData.m_xMin = xMin;
    clusterData.m_xMax = xMax;
    clusterData.m_wireMin = wireMin;
    clusterData.m_wireMax = wireMax;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// GroupClustersByOverlap
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::GroupClustersByOverlap(
    const std::vector<ClusterData> &allClusters, std::vector<ClusterGroup> &groups) const
{
    if (allClusters.empty())
        return;

    // Collect all cluster X-range endpoints and merge overlapping intervals.
    // A valid cut boundary can only exist in a gap where no cluster is present
    // — i.e. no cluster has xMin < cut < xMax. Merging all overlapping X-ranges
    // gives exactly the set of regions that must stay together.
    std::vector<std::pair<float, float>> intervals;
    intervals.reserve(allClusters.size());

    for (const ClusterData &cd : allClusters)
        intervals.emplace_back(cd.m_xMin, cd.m_xMax);

    std::sort(intervals.begin(), intervals.end());

    std::vector<std::pair<float, float>> merged;
    for (const auto &iv : intervals)
    {
        if (!merged.empty() && iv.first <= merged.back().second + m_xOverlapThreshold)
            merged.back().second = std::max(merged.back().second, iv.second);
        else
            merged.push_back(iv);
    }

    // Assign each cluster to the merged interval that contains its X range.
    // By construction every cluster falls entirely within exactly one merged interval.
    std::map<std::size_t, ClusterGroup> groupMap;

    for (const ClusterData &cd : allClusters)
    {
        for (std::size_t m = 0; m < merged.size(); ++m)
        {
            if (cd.m_xMin >= merged[m].first && cd.m_xMax <= merged[m].second)
            {
                ClusterGroup &group(groupMap[m]);
                group.m_clusters.push_back(cd);
                group.m_xOverlapMin = merged[m].first;
                group.m_xOverlapMax = merged[m].second;
                break;
            }
        }
    }

    for (auto &kv : groupMap)
        groups.push_back(std::move(kv.second));
}

//------------------------------------------------------------------------------------------------------------------------------------------
// ComputeProjectionMatrix
//------------------------------------------------------------------------------------------------------------------------------------------

GaussianSplattingAlgorithm::ProjectionMatrix GaussianSplattingAlgorithm::ComputeProjectionMatrix(
    const HitType hitType) const
{
    // P_k maps (x, y, z) -> (x, wire_coord).
    // wire_coord is a linear function of (y, z) only; x is preserved.
    // Probe LArGeometryHelper with unit vectors to extract coefficients.
    const CartesianVector origin(0.f, 0.f, 0.f);
    const CartesianVector ey(0.f, 1.f, 0.f);
    const CartesianVector ez(0.f, 0.f, 1.f);

    const CartesianVector p0(LArGeometryHelper::ProjectPosition(this->GetPandora(), origin, hitType));
    const CartesianVector py(LArGeometryHelper::ProjectPosition(this->GetPandora(), ey, hitType));
    const CartesianVector pz(LArGeometryHelper::ProjectPosition(this->GetPandora(), ez, hitType));

    // Partial derivatives of wire_coord w.r.t y and z
    const double dWire_dy(static_cast<double>(py.GetZ() - p0.GetZ()));
    const double dWire_dz(static_cast<double>(pz.GetZ() - p0.GetZ()));

    ProjectionMatrix proj;
    // Row 0: x is preserved as-is
    // Row 1: wire coordinate = offset + dWire_dy*y + dWire_dz*z
    //        (offset absorbed when computing δ = hit - projected_mean)
    proj.m_P << 1.0, 0.0, 0.0,
                0.0, dWire_dy, dWire_dz;
    return proj;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// InitialiseGaussians
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::InitialiseGaussians(
    const ClusterGroup &group, const std::map<HitType, ProjectionMatrix> &projections, std::vector<Gaussian3D> &gaussians) const
{
    for (const ClusterData &cd : group.m_clusters)
    {
        const ProjectionMatrix &P(projections.at(cd.m_hitType));
        const double dW_dy(P.m_P(1, 1));
        const double dW_dz(P.m_P(1, 2));

        for (const Hit2D &hit : cd.m_hits)
        {
            const double xHit(hit.m_position(0));
            const double wireHit(hit.m_position(1));

            // Place mean at the hit's known (x, wire) position.
            // The ambiguous 3D coordinate is unconstrained by this view alone;
            // initialise it to 0 (detector midpoint) and solve the wire equation
            // for the constrained direction.
            double initY(0.0), initZ(0.0);

            if (std::abs(dW_dz) > std::abs(dW_dy))
            {
                initY = 0.0;
                initZ = (std::abs(dW_dz) > 1e-9) ? (wireHit / dW_dz) : 0.0;
            }
            else
            {
                initZ = 0.0;
                initY = (std::abs(dW_dy) > 1e-9) ? (wireHit / dW_dy) : 0.0;
            }

            // Sigma in x: match the hit's drift sampling width.
            // Sigma in the wire-constrained direction: back-project hit wire width.
            //   projected wire sigma = |dW_dq| * sigma_q  for the dominant direction q.
            // Sigma in the ambiguous direction: large initial value (will be shaped by
            // constraints from other views during optimisation).
            double sigmaQ(static_cast<double>(m_initialSigmaXY));  // ambiguous dimension
            double sigmaY(sigmaQ), sigmaZ(sigmaQ);

            if (std::abs(dW_dz) > std::abs(dW_dy))
                sigmaZ = (std::abs(dW_dz) > 1e-9) ? (hit.m_sigmaWire / std::abs(dW_dz)) : sigmaQ;
            else
                sigmaY = (std::abs(dW_dy) > 1e-9) ? (hit.m_sigmaWire / std::abs(dW_dy)) : sigmaQ;

            Gaussian3D g;
            g.m_mean     = Eigen::Vector3d(xHit, initY, initZ);
            g.m_logSigma = Eigen::Vector3d(std::log(hit.m_sigmaX), std::log(sigmaY), std::log(sigmaZ));
            g.m_opacity  = 1.0;

            gaussians.push_back(g);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// RunGradientDescent
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::RunGradientDescent(
    const ClusterGroup &group, const std::map<HitType, ProjectionMatrix> &projections, std::vector<Gaussian3D> &gaussians) const
{
    const std::size_t nGaussians(gaussians.size());

    // Initialise Adam states: 7 parameters per Gaussian
    // (3 mean + 3 logSigma + 1 opacity)
    std::vector<AdamState> adamStates(nGaussians);

    for (AdamState &state : adamStates)
    {
        state.m_firstMoment = Eigen::VectorXd::Zero(7);
        state.m_secondMoment = Eigen::VectorXd::Zero(7);
        state.m_step = 0;
    }

    std::vector<Eigen::Vector3d> gradMean(nGaussians);
    std::vector<Eigen::Vector3d> gradLogSigma(nGaussians);
    std::vector<double> gradOpacity(nGaussians);

    for (unsigned int step = 0; step < m_nOptimisationSteps; ++step)
    {
        this->ComputeLossAndGradients(group, projections, gaussians, gradMean, gradLogSigma, gradOpacity);
        this->ApplyAdamStep(gaussians, adamStates, gradMean, gradLogSigma, gradOpacity);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// ComputeLossAndGradients
//------------------------------------------------------------------------------------------------------------------------------------------

double GaussianSplattingAlgorithm::ComputeLossAndGradients(const ClusterGroup &group,
    const std::map<HitType, ProjectionMatrix> &projections, const std::vector<Gaussian3D> &gaussians,
    std::vector<Eigen::Vector3d> &gradMean, std::vector<Eigen::Vector3d> &gradLogSigma, std::vector<double> &gradOpacity) const
{
    const std::size_t nGaussians(gaussians.size());
    double totalLoss(0.0);

    // Zero all gradients
    for (std::size_t g = 0; g < nGaussians; ++g)
    {
        gradMean[g].setZero();
        gradLogSigma[g].setZero();
        gradOpacity[g] = 0.0;
    }

    // Iterate over every view's clusters and their hits
    for (const ClusterData &cd : group.m_clusters)
    {
        const ProjectionMatrix &proj(projections.at(cd.m_hitType));
        const Eigen::Matrix<double, 2, 3> &P(proj.m_P);

        // P row 1 coefficients for y and z
        const double dW_dy(P(1, 1));
        const double dW_dz(P(1, 2));

        for (const Hit2D &hit : cd.m_hits)
        {
            const double hx(hit.m_position(0));
            const double hz(hit.m_position(1)); // wire coord
            const double sigHitX(hit.m_sigmaX);
            const double sigHitWire(hit.m_sigmaWire);

            for (std::size_t g = 0; g < nGaussians; ++g)
            {
                const Gaussian3D &gauss(gaussians[g]);

                // Projected mean: (μ_x, dW_dy*μ_y + dW_dz*μ_z)
                const double muX(gauss.m_mean(0));
                const double muWire(dW_dy * gauss.m_mean(1) + dW_dz * gauss.m_mean(2));

                // Gaussian sigma in each projected direction:
                // sX² = σ_x² + σ_hit_x²
                // sWire² = (dW_dy*σ_y)² + (dW_dz*σ_z)² + σ_hit_wire²
                const Eigen::Vector3d sigma = gauss.m_logSigma.array().exp();
                const double sigX2(sigma(0) * sigma(0) + sigHitX * sigHitX);
                const double sigWire2(dW_dy * dW_dy * sigma(1) * sigma(1) + dW_dz * dW_dz * sigma(2) * sigma(2) +
                    sigHitWire * sigHitWire);

                // Residuals
                const double dX(hx - muX);
                const double dW(hz - muWire);

                // Gaussian density: N(h; μ_k, S_k) = (1/(2π√(sX²·sWire²))) exp(-0.5*(dX²/sX² + dW²/sWire²))
                const double exponent(-0.5 * (dX * dX / sigX2 + dW * dW / sigWire2));

                // Guard against underflow
                if (exponent < -50.0)
                    continue;

                const double normFactor(1.0 / (2.0 * M_PI * std::sqrt(sigX2 * sigWire2)));
                const double density(normFactor * std::exp(exponent));
                const double f(gauss.m_opacity * density); // α_g * N

                totalLoss -= f;

                // -------------------------------------------------------
                // Gradients of f = α_g * N w.r.t. Gaussian parameters
                // ∂f/∂μ_x    = f * (dX / sX²)
                // ∂f/∂μ_y    = f * (dW / sWire²) * dW_dy
                // ∂f/∂μ_z    = f * (dW / sWire²) * dW_dz
                // ∂f/∂l_x    = f * (dX²/sX² - 1) * σ_x²/sX²   [chain through exp(2*l_x)]
                // ∂f/∂l_y    = f * (dW²/sWire² - 1) * (dW_dy*σ_y)²/sWire²
                // ∂f/∂l_z    = f * (dW²/sWire² - 1) * (dW_dz*σ_z)²/sWire²
                // ∂f/∂α_g    = N (= f / α_g)
                // Loss = -Σ f  =>  ∂L/∂p = -∂f/∂p
                // -------------------------------------------------------

                const double dX_over_sX2(dX / sigX2);
                const double dW_over_sWire2(dW / sigWire2);

                gradMean[g](0) -= f * dX_over_sX2;
                gradMean[g](1) -= f * dW_over_sWire2 * dW_dy;
                gradMean[g](2) -= f * dW_over_sWire2 * dW_dz;

                const double chiX2(dX * dX / sigX2);
                const double chiWire2(dW * dW / sigWire2);
                const double sigX2_g(sigma(0) * sigma(0));
                const double sigY2_g(sigma(1) * sigma(1));
                const double sigZ2_g(sigma(2) * sigma(2));

                gradLogSigma[g](0) -= f * (chiX2 - 1.0) * sigX2_g / sigX2;
                gradLogSigma[g](1) -= f * (chiWire2 - 1.0) * (dW_dy * dW_dy * sigY2_g) / sigWire2;
                gradLogSigma[g](2) -= f * (chiWire2 - 1.0) * (dW_dz * dW_dz * sigZ2_g) / sigWire2;

                gradOpacity[g] -= density;
            }
        }
    }

    return totalLoss;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// ApplyAdamStep
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::ApplyAdamStep(std::vector<Gaussian3D> &gaussians, std::vector<AdamState> &adamStates,
    const std::vector<Eigen::Vector3d> &gradMean, const std::vector<Eigen::Vector3d> &gradLogSigma,
    const std::vector<double> &gradOpacity) const
{
    const double lr(static_cast<double>(m_learningRate));
    const double b1(static_cast<double>(m_adamBeta1));
    const double b2(static_cast<double>(m_adamBeta2));
    const double eps(static_cast<double>(m_adamEpsilon));

    for (std::size_t g = 0; g < gaussians.size(); ++g)
    {
        AdamState &state(adamStates[g]);
        ++(state.m_step);
        const int t(state.m_step);

        // Pack gradients into a single 7-vector: [mean(3), logSigma(3), opacity(1)]
        Eigen::VectorXd grad(7);
        grad << gradMean[g](0), gradMean[g](1), gradMean[g](2),
                gradLogSigma[g](0), gradLogSigma[g](1), gradLogSigma[g](2),
                gradOpacity[g];

        // Bias-corrected moment estimates
        state.m_firstMoment  = b1 * state.m_firstMoment + (1.0 - b1) * grad;
        state.m_secondMoment = b2 * state.m_secondMoment + (1.0 - b2) * grad.array().square().matrix();

        const double biasCorr1(1.0 - std::pow(b1, t));
        const double biasCorr2(1.0 - std::pow(b2, t));

        const Eigen::VectorXd mHat(state.m_firstMoment / biasCorr1);
        const Eigen::VectorXd vHat(state.m_secondMoment / biasCorr2);

        // Parameter update: θ ← θ - lr * m̂ / (√v̂ + ε)
        const Eigen::VectorXd update(lr * mHat.array() / (vHat.array().sqrt() + eps));

        Gaussian3D &gauss(gaussians[g]);
        gauss.m_mean(0) -= update(0);
        gauss.m_mean(1) -= update(1);
        gauss.m_mean(2) -= update(2);
        gauss.m_logSigma(0) -= update(3);
        gauss.m_logSigma(1) -= update(4);
        gauss.m_logSigma(2) -= update(5);

        // Opacity update: clamp to (pruneThreshold/2, 1)
        gauss.m_opacity = std::max(static_cast<double>(m_opacityPruneThreshold) * 0.5,
            std::min(1.0, gauss.m_opacity - update(6)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// PruneGaussians
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::PruneGaussians(std::vector<Gaussian3D> &gaussians) const
{
    gaussians.erase(
        std::remove_if(gaussians.begin(), gaussians.end(),
            [this](const Gaussian3D &g) { return g.m_opacity < static_cast<double>(m_opacityPruneThreshold); }),
        gaussians.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------
// VisualiseGaussians
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::VisualiseGaussians(
    const std::vector<Gaussian3D> &gaussians, const std::map<pandora::HitType, ProjectionMatrix> &projections) const
{
#ifdef MONITORING
    // View types and their display names
    const std::vector<std::pair<HitType, std::string>> views = {
        {TPC_VIEW_U, "U"}, {TPC_VIEW_V, "V"}, {TPC_VIEW_W, "W"}};

    // Cycle colours across Gaussians so they can be distinguished
    const std::vector<Color> colors = {RED, GREEN, BLUE, ORANGE, MAGENTA, CYAN, YELLOW, WHITE};

    for (std::size_t gi = 0; gi < gaussians.size(); ++gi)
    {
        const Gaussian3D &g(gaussians[gi]);
        const Color color(colors[gi % colors.size()]);
        const Eigen::Vector3d sigma(g.m_logSigma.array().exp());

        for (const auto &[hitType, viewName] : views)
        {
            const auto projIt(projections.find(hitType));
            if (projIt == projections.end())
                continue;
            const ProjectionMatrix &P(projIt->second);

            // Project mean: (x, wire) = P.m_P * mu
            const Eigen::Vector2d projMu(P.m_P * g.m_mean);
            const float cx(static_cast<float>(projMu(0)));
            const float cw(static_cast<float>(projMu(1)));

            // Project covariance: S = P.m_P * diag(sigma^2) * P.m_P^T  (2x2)
            const Eigen::Matrix3d Sigma3(Eigen::Vector3d(sigma.array().square()).asDiagonal());
            const Eigen::Matrix2d S(P.m_P * Sigma3 * P.m_P.transpose());

            // Eigen-decompose S to get major/minor axes
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(S);
            // Eigenvalues in ascending order; take sqrt for sigma lengths
            const double sigma0(std::sqrt(std::max(0.0, eig.eigenvalues()(0))));
            const double sigma1(std::sqrt(std::max(0.0, eig.eigenvalues()(1))));
            const Eigen::Vector2d axis0(eig.eigenvectors().col(0)); // minor axis
            const Eigen::Vector2d axis1(eig.eigenvectors().col(1)); // major axis

            // Draw major axis line: centre ± 1σ along axis1
            // In 2D view: x is drift, z is wire coordinate
            const CartesianVector majStart(cx - static_cast<float>(sigma1 * axis1(0)), 0.f,
                cw - static_cast<float>(sigma1 * axis1(1)));
            const CartesianVector majEnd(cx + static_cast<float>(sigma1 * axis1(0)), 0.f,
                cw + static_cast<float>(sigma1 * axis1(1)));

            // Draw minor axis line: centre ± 1σ along axis0
            const CartesianVector minStart(cx - static_cast<float>(sigma0 * axis0(0)), 0.f,
                cw - static_cast<float>(sigma0 * axis0(1)));
            const CartesianVector minEnd(cx + static_cast<float>(sigma0 * axis0(0)), 0.f,
                cw + static_cast<float>(sigma0 * axis0(1)));

            const std::string labelMaj("G" + std::to_string(gi) + "_" + viewName + "_maj");
            const std::string labelMin("G" + std::to_string(gi) + "_" + viewName + "_min");

            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &majStart, &majEnd, labelMaj, color, 2, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &minStart, &minEnd, labelMin, color, 1, 2));
        }
    }
#else
    (void)gaussians;
    (void)projections;
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------
// ComputeResponsibility
//------------------------------------------------------------------------------------------------------------------------------------------

double GaussianSplattingAlgorithm::ComputeResponsibility(
    const Gaussian3D &g, const ClusterData &cd, const ProjectionMatrix &projection) const
{
    // Compute r(g, c_k) = P(Gaussian g's projected mass falls in cluster c_k's bounding box)
    //
    // For the projected 2D Gaussian:
    //   μ_k = (μ_x, dW_dy*μ_y + dW_dz*μ_z)
    //   sX²  = σ_x²          (x projected trivially, no hit uncertainty here)
    //   sWire² = (dW_dy*σ_y)² + (dW_dz*σ_z)²
    //
    // r(g, c_k) = Φ_x(xMax, xMin) * Φ_w(wireMax, wireMin)
    // where Φ_x(a,b) = Φ((a-μ_x)/sX) - Φ((b-μ_x)/sX)  (Gaussian CDF difference)

    const Eigen::Matrix<double, 2, 3> &P(projection.m_P);
    const double dW_dy(P(1, 1));
    const double dW_dz(P(1, 2));

    const Eigen::Vector3d sigma(g.m_logSigma.array().exp());

    const double muX(g.m_mean(0));
    const double muWire(dW_dy * g.m_mean(1) + dW_dz * g.m_mean(2));

    const double sX(sigma(0));
    const double sWire(std::sqrt(dW_dy * dW_dy * sigma(1) * sigma(1) + dW_dz * dW_dz * sigma(2) * sigma(2)));

    if (sX < 1e-9 || sWire < 1e-9)
        return 0.0;

    // Gaussian CDF using erfc: Φ(z) = 0.5 * erfc(-z / sqrt(2))
    auto gaussianCDF = [](double x, double mu, double sigma_val) -> double
    {
        return 0.5 * std::erfc(-(x - mu) / (sigma_val * M_SQRT2));
    };

    const double cdfXHigh(gaussianCDF(static_cast<double>(cd.m_xMax), muX, sX));
    const double cdfXLow(gaussianCDF(static_cast<double>(cd.m_xMin), muX, sX));
    const double cdfWHigh(gaussianCDF(static_cast<double>(cd.m_wireMax), muWire, sWire));
    const double cdfWLow(gaussianCDF(static_cast<double>(cd.m_wireMin), muWire, sWire));

    const double fracX(std::max(0.0, cdfXHigh - cdfXLow));
    const double fracWire(std::max(0.0, cdfWHigh - cdfWLow));

    return fracX * fracWire;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// ComputeCoverageScores
//------------------------------------------------------------------------------------------------------------------------------------------

void GaussianSplattingAlgorithm::ComputeCoverageScores(const ClusterGroup &group,
    const std::map<HitType, ProjectionMatrix> &projections, const std::vector<Gaussian3D> &gaussians,
    CoverageScoreMap &scoreMap) const
{
    const std::size_t nClusters(group.m_clusters.size());
    const std::size_t nGaussians(gaussians.size());

    // Pre-compute r(g, c_k) for all (Gaussian, cluster) pairs
    // responsibilities[g][k] = r(g, c_k)
    std::vector<std::vector<double>> responsibilities(nGaussians, std::vector<double>(nClusters, 0.0));

    for (std::size_t g = 0; g < nGaussians; ++g)
    {
        for (std::size_t k = 0; k < nClusters; ++k)
        {
            const ClusterData &cd(group.m_clusters[k]);
            const ProjectionMatrix &proj(projections.at(cd.m_hitType));
            responsibilities[g][k] = this->ComputeResponsibility(gaussians[g], cd, proj);
        }
    }

    // Pre-compute per-Gaussian weight: α_g * V_g
    std::vector<double> gaussianWeight(nGaussians);

    for (std::size_t g = 0; g < nGaussians; ++g)
        gaussianWeight[g] = gaussians[g].m_opacity * gaussians[g].GetVolume();

    // Compute S(c_i, c_j) and S_self(c_i) for all pairs of clusters from
    // different views only (cross-view affinities are what link 2D to 3D)
    for (std::size_t i = 0; i < nClusters; ++i)
    {
        const ClusterData &cdi(group.m_clusters[i]);

        // S_self(c_i) = Σ_g r(g,c_i)^2 * α_g * V_g
        double selfI(0.0);

        for (std::size_t g = 0; g < nGaussians; ++g)
            selfI += responsibilities[g][i] * responsibilities[g][i] * gaussianWeight[g];

        for (std::size_t j = i + 1; j < nClusters; ++j)
        {
            const ClusterData &cdj(group.m_clusters[j]);

            // Only score cross-view pairs
            if (cdi.m_hitType == cdj.m_hitType)
                continue;

            // S_self(c_j)
            double selfJ(0.0);

            for (std::size_t g = 0; g < nGaussians; ++g)
                selfJ += responsibilities[g][j] * responsibilities[g][j] * gaussianWeight[g];

            // S(c_i, c_j) = Σ_g r(g,c_i) * r(g,c_j) * α_g * V_g
            double shared(0.0);

            for (std::size_t g = 0; g < nGaussians; ++g)
                shared += responsibilities[g][i] * responsibilities[g][j] * gaussianWeight[g];

            // F(c_i -> c_j) = S(c_i,c_j) / S_self(c_i)
            // F(c_j -> c_i) = S(c_i,c_j) / S_self(c_j)
            const double fIJ((selfI > 1e-30) ? (shared / selfI) : 0.0);
            const double fJI((selfJ > 1e-30) ? (shared / selfJ) : 0.0);

            CoverageScores scores;
            scores.m_forwardCoverage = static_cast<float>(std::min(fIJ, 1.0));
            scores.m_reverseCoverage = static_cast<float>(std::min(fJI, 1.0));

            scoreMap[{cdi.m_pCluster, cdj.m_pCluster}] = scores;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// ReadSettings
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GaussianSplattingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "XOverlapThreshold", m_xOverlapThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NumGradientSteps", m_nOptimisationSteps));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NOptimisationSteps", m_nOptimisationSteps));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "LearningRate", m_learningRate));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AdamBeta1", m_adamBeta1));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AdamBeta2", m_adamBeta2));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AdamEpsilon", m_adamEpsilon));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InitialSigmaXY", m_initialSigmaXY));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InitialSigmaWire", m_initialSigmaWire));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OpacityThreshold", m_opacityPruneThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OpacityPruneThreshold", m_opacityPruneThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchThreshold", m_matchThreshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
