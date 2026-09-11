/**
 *  @file   larpandoracontent/LArSolutions/SolutionSimpleClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the SolutionSimpleClusterSplittingAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SOLUTION_SIMPLE_CLUSTER_SPLITTING_H
#define LAR_SOLUTION_SIMPLE_CLUSTER_SPLITTING_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SolutionSimpleClusterSplittingAlgorithm class
 */
class SolutionSimpleClusterSplittingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SolutionSimpleClusterSplittingAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetMCMuon(const pandora::MCParticleList *const pMCParticleList, const pandora::MCParticle *&pMCMuon) const;
    
    void SplitCluster(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition) const;    

    bool GetPathwayDirections(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition, pandora::CartesianVector &splitAxis1, pandora::CartesianVector &splitAxis2) const;
    
    void SplitCaloHits(const pandora::Cluster *const pCluster, const pandora::CartesianVector &splitPosition, const pandora::CartesianVector &splitAxis1,
        const pandora::CartesianVector &splitAxis2, pandora::CaloHitList &caloHitList1, pandora::CaloHitList &caloHitList2) const;
    
    std::string m_clusterListName;
    std::string m_mcParticleListName;
    int m_nAngularBins;
    float m_radiusForDirEstimate;
    float m_maxDistToSplitPos;
};

} // namespace lar_content

#endif // #ifndef LAR_SOLUTION_SIMPLE_CLUSTER_SPLITTING_H

