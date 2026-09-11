/**
 *  @file   larpandoracontent/LArSolution/SolutionLightClusteringAlgorithm.h
 *
 *  @brief  Header file for the SolutionLightClusteringAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SOLUTION_LIGHT_CLUSTERING_H
#define LAR_SOLUTION_LIGHT_CLUSTERING_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

namespace lar_content
{

/**
 *  @brief  SolutionLightClusteringAlgorithm class
 */
class SolutionLightClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SolutionLightClusteringAlgorithm();

private:
    typedef std::vector<const LArOpHit *> OpHitVector;
    
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void PrepareOpHitList(const pandora::CaloHitList *const pOpHitList, OpHitVector &opHitVector) const;

    bool GetSeed(const OpHitVector &opHitVector, const OpHitVector &consideredSeeds, const LArOpHit *& pOpSeed) const;

    void GetClusterHits(const LArOpHit *const pOpSeed, OpHitVector &opHitVector, pandora::CaloHitList &opClusterHits) const;

    bool IsClose(const LArOpHit *const pEventOpHit, const LArOpHit *const pClusterOpHit) const;
    
    std::string m_opticalHitListName;
    std::string m_outputOpClusterListName;
    float m_opHitPEThreshold;
    float m_maxDeltaY;
    float m_maxDeltaZ;
    float m_maxDeltaT;    
};

} // namespace lar_content

#endif // #ifndef LAR_SOLUTION_LIGHT_CLUSTERING_H

