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
    typedef std::map<const pandora::Cluster *const, std::map<const pandora::Cluster *const, float>> ClusterRelationshipMatrix;
    typedef std::map<const pandora::MCParticle *const, int> MCCountsMap;
    typedef std::map<const pandora::MCParticle *const, pandora::CaloHitList> MCCaloHitsMap;
    typedef std::map<const pandora::Cluster *const, std::vector<const pandora::Cluster *>> AdjacencyLists;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const pandora::MCParticle *GetMainMC(const pandora::CaloHit *const pCaloHit) const;

    const pandora::MCParticle *FoldMCTo(const pandora::MCParticle *const pMC) const;

    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendentElectrons) const;

    const pandora::MCParticle *FoldPotentialDeltaRayTo(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pMC) const;

    pandora::StatusCode PerformIntraClusterMerges() const;

    pandora::StatusCode CalcTargetClusterRelationships(ClusterRelationshipMatrix &clusterRelationshipMatrix) const;

    pandora::StatusCode GetClusters(pandora::ClusterList &clusterList) const;

    float CalcIntraViewClusterRelationship(
        const MCCountsMap &mcCountsI, const int totalCountsI, const MCCountsMap &mcCountsJ, const int totalCountsJ) const;

    float CalcInterViewClusterRelationship(
        const MCCaloHitsMap &mcCaloHitsI, const int totalCountsI, const MCCaloHitsMap &mcCaloHitsJ) const;

    pandora::StatusCode GetClusterMCSummary(
        const pandora::Cluster *const pCluster, MCCountsMap &mcCounts, int &totalCounts, MCCaloHitsMap &mcCaloHits) const;

    int CalcMaxXOverlapHitMatches(const pandora::CaloHitList &caloHitsI, const pandora::CaloHitList &caloHitsJ) const;

    bool FindAugmentingPath(
        const int idxI, const std::vector<std::vector<int>> &adjIToJ, std::vector<int> &matchesJ, std::vector<bool> &visitedJ) const;

    pandora::StatusCode PopulateAdjacencyLists(
        const ClusterRelationshipMatrix &relMat, AdjacencyLists &coreClusterAdjLists, AdjacencyLists &accClusterAdjLists) const;

    pandora::StatusCode CalculateConnectedGroups(
        const AdjacencyLists &clusterAdjLists, std::vector<pandora::ClusterSet> &clusterGroups) const;

    pandora::StatusCode AddClustersToGroups(
        const ClusterRelationshipMatrix &clusterRelMat,
        AdjacencyLists &ungroupedClusterAdjLists,
        std::vector<pandora::ClusterSet> &clusterGroups) const;

    bool IsSingletonPartition(const std::vector<pandora::ClusterSet> &clusterGroups) const;

    pandora::StatusCode MergeGroups(const std::vector<pandora::ClusterSet> &clusterGroups, const std::string &listNmae) const;

    /* Start shared mutable members */

    mutable std::map<const pandora::MCParticle *const, const pandora::MCParticle *const> m_mcFoldTo; // caching

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
    float m_scoreThreshold;
    unsigned int m_accessoryClustersMaxHits;

    /* End configurable via xml members */
};

} // namespace lar_dl_content

#endif // LAR_DL_MATCHING_ALGORITHM
