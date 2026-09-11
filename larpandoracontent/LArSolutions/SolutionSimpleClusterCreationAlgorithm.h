/**
 *  @file   larpandoracontent/LArSolutions/SolutionSimpleClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the SolutionSimpleClusterCreationAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SOLUTION_SIMPLE_CLUSTER_CREATION_H
#define LAR_SOLUTION_SIMPLE_CLUSTER_CREATION_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SolutionSimpleClusterCreationAlgorithm class
 */
class SolutionSimpleClusterCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SolutionSimpleClusterCreationAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetClusterHits(const pandora::CaloHit *const pSeedHit, const pandora::CaloHitVector &eventHits, pandora::CaloHitList &clusterHits) const;    

    std::string m_caloHitListName;
    std::string m_outputClusterListName;
    unsigned int m_maxNCaloHits;
    float m_maxSeparation;
};

} // namespace lar_content

#endif // #ifndef LAR_SOLUTION_SIMPLE_CLUSTER_CREATION_H

