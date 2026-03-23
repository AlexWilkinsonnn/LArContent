/**
 *  @file   
 *
 *  @brief  
 *
 *  $Log: $
 */
#ifndef LAR_DL_MATCHING_ALGORITHM
#define LAR_DL_MATCHING_ALGORITHM 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

namespace lar_dl_content
{
/**
 *  @brief  DLMatchingAlgorithm class
 */
class DLMatchingAlgorithm : public pandora::Algorithm
{
private:

public:
    /**
     *  @brief Default constructor
     */
    DLMatchingAlgorithm();

    /**
     *  @brief Default destructor
     */
    virtual ~DLMatchingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const pandora::MCParticle *GetMainMC(const pandora::CaloHit *const pCaloHit) const;

    const pandora::MCParticle *FoldMCTo(const pandora::MCParticle *const pMC) const;

    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendentElectrons) const;

    const pandora::MCParticle *FoldPotentialDeltaRayTo(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pMC) const;

    pandora::StatusCode CalcTargetClusterRelationships(
        std::map<const pandora::Cluster *const, std::map<const pandora::Cluster *const, float>> &clusterRelationshipMatrix) const;

    pandora::StatusCode GetClusters(pandora::ClusterList &clusterList) const;

    float CalcIntraViewClusterRelationship(
        const std::map<const pandora::MCParticle *const, int> &mcCountsI, const int totalCountsI, 
        const std::map<const pandora::MCParticle *const, int> &mcCountsJ, const int totalCountsJ) const;

    float CalcInterViewClusterRelationship(
        const std::map<const pandora::MCParticle *const, pandora::CaloHitList> &mcCaloHitsI, const int totalCountsI,
        const std::map<const pandora::MCParticle *const, pandora::CaloHitList> &mcCaloHitsJ) const;

    pandora::StatusCode GetClusterMCSummary(
        const pandora::Cluster *const pCluster,
        std::map<const pandora::MCParticle *const, int> &mcCounts, int &totalCounts,
        std::map<const pandora::MCParticle *const, pandora::CaloHitList> &mcCaloHits) const;

    int CalcMaxXOverlapHitMatches(const pandora::CaloHitList &caloHitsI, const pandora::CaloHitList &caloHitsJ) const;

    bool FindAugmentingPath(
        const int idxI, const std::vector<std::vector<int>> &adjIToJ, std::vector<int> &matchesJ, std::vector<bool> &visitedJ) const;

    /* Start shared mutable members */

    mutable std::map<const pandora::MCParticle *const, const pandora::MCParticle *const> m_mcFoldTo;

    /* End shared mutable members */

    /* Start hardcoded members */

    std::map<pandora::HitType, float> m_deltaRayLengthThresholdSquared;
    float m_deltaRayParentWeightThreshold;

    /* End hardcoded members */

    /* Start configurable via xml members */

    std::string m_inputClusterListNameU;
    std::string m_inputClusterListNameV;
    std::string m_inputClusterListNameW;
    float m_xOverlapMatchingLeeway;

    /* End configurable via xml members */
};

} // namespace lar_dl_content

#endif // LAR_DL_MATCHING_ALGORITHM
