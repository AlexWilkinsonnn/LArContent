/**
 *  @file   larpandoracontent//SolutionExampleAlgorithm.cc
 *
 *  @brief  Implementation of the SolutionExampleAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArSolutions/SolutionExampleAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SolutionExampleAlgorithm::SolutionExampleAlgorithm() :
    m_favouriteNumber(21)
{
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionExampleAlgorithm::Run()
{
    std::cout << "The best number is: " << m_favouriteNumber << std::endl;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    std::cout << "pMCParticleList->size(): " << pMCParticleList->size() << std::endl;
    
    MCParticleVector mcParticleVec(pMCParticleList->begin(), pMCParticleList->end());
    std::sort(mcParticleVec.begin(), mcParticleVec.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        std::cout << "PDG: " << pMCParticle->GetParticleId() << std::endl;
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionExampleAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FavouriteNumber", m_favouriteNumber));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

