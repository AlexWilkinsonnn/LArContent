/**
 *  @file   larpandoracontent/LArSolutions/SolutionLightAlgorithm.h
 *
 *  @brief  Header file for the SolutionLightAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SOLUTION_LIGHT_H
#define LAR_SOLUTION_LIGHT_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"

namespace lar_content
{

/**
 *  @brief  SolutionLightAlgorithm class
 */
class SolutionLightAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SolutionLightAlgorithm();

    ~SolutionLightAlgorithm() override;

private:
    pandora::StatusCode Run() override;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    bool IsMichelElectron(const LArMCParticle *const pMC) const;

    std::string m_opticalHitListName;
    std::string m_fileName;
    std::string m_treeName;
};

} // namespace lar_content

#endif // #ifndef LAR_SOLUTION_LIGHT_H

