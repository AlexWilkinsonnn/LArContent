/**
 *  @file   larpandoracontent/LArSolutions/SolutionChargeLightMatchingAlgorithm.h
 *
 *  @brief  Header file for the SolutionChargeLightMatchingAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SOLUTION_CHARGE_LIGHT_MATCHING_H
#define LAR_SOLUTION_CHARGE_LIGHT_MATCHING_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SolutionChargeLightMatchingAlgorithm class
 */
class SolutionChargeLightMatchingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SolutionChargeLightMatchingAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetClosestCaloHitToPosition(const pandora::CartesianVector &position, const pandora::CaloHitList &caloHitList, const pandora::CaloHit *&pClosestHit);
    
    std::string m_opticalClusterListName;
    std::string m_pfoListName;
};

} // namespace lar_content

#endif // #ifndef LAR_SOLUTION_CHARGE_LIGHT_MATCHING_H

