/**
 *  @file   larpandoracontent/LArSolutions/SolutionLightAlgorithm.cc
 *
 *  @brief  Implementation of the SolutionLightAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArSolutions/SolutionLightAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

SolutionLightAlgorithm::SolutionLightAlgorithm() :
    m_fileName("LightValidation.root"),
    m_treeName("events")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SolutionLightAlgorithm::~SolutionLightAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "RECREATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionLightAlgorithm::Run()
{
    const CaloHitList *pOpHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_opticalHitListName, pOpHitList));

    std::vector<float> opHitPEs, opHitPeakTimes, opHitStartTimes, opHitWidths, opHitPositionsX, opHitPositionsY, opHitPositionsZ;
    std::vector<int> opHitChannels;

    std::cout << "Got " << m_opticalHitListName << " list with " << pOpHitList->size() << " hits" << std::endl;
    // int linesPrinted{0};
    for (const CaloHit *const pCaloHit : *pOpHitList)
    {
        const LArOpHit *const pOpHit(dynamic_cast<const LArOpHit *>(pCaloHit));
        if (!pOpHit)
            return STATUS_CODE_INVALID_PARAMETER;

        // if (++linesPrinted < 5)
        // {
        //   std::cout << "-- Op Hit " << linesPrinted << "\n"
        //             << "---- channel = " << pOpHit->GetChannel() << " at position " << pOpHit->GetPositionVector() << "\n"
        //             << "---- peak time = " << pOpHit->GetTime() << " us\n"
        //             << "---- start time = " << pOpHit->GetStartTime() << " us\n"
        //             << "---- end time = " << pOpHit->GetStartTime() + pOpHit->GetWidth() << " us\n"
        //             << "---- PE = " << pOpHit->GetInputEnergy() << std::endl;
        // }
        // else if (linesPrinted == 5)
        // {
        //     std::cout << "-- ..." << std::endl;
        // }

        opHitPEs.push_back(pOpHit->GetInputEnergy());
        opHitPeakTimes.push_back(pOpHit->GetTime());
        opHitStartTimes.push_back(pOpHit->GetStartTime());
        opHitWidths.push_back(pOpHit->GetWidth());
        opHitPositionsX.push_back(pOpHit->GetPositionVector().GetX());
        opHitPositionsY.push_back(pOpHit->GetPositionVector().GetY());
        opHitPositionsZ.push_back(pOpHit->GetPositionVector().GetZ());
        opHitChannels.push_back(static_cast<int>(pOpHit->GetChannel()));
    }

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    std::vector<float> mcParticleTimes, mcParticleEnergies;
    std::vector<int> mcParticlePdgs, mcParticleIsMichel;

    std::vector<float> mcParticleTrajX;
    std::vector<float> mcParticleTrajY;
    std::vector<float> mcParticleTrajZ;
    std::vector<float> mcParticleTrajT;

    for (const MCParticle *const pMC : *pMCParticleList)
    {
        const LArMCParticle *const pLArMC(dynamic_cast<const LArMCParticle *>(pMC));
        if (!pLArMC)
            return STATUS_CODE_INVALID_PARAMETER;
        mcParticlePdgs.push_back(pLArMC->GetParticleId());
        mcParticleTimes.push_back(pLArMC->GetT0());
        mcParticleEnergies.push_back(pLArMC->GetEnergy());
        mcParticleIsMichel.push_back(this->IsMichelElectron(pLArMC) ? 1 : 0);

        CartesianPointVector trajPoints(pLArMC->GetTrajPoints());
        for (const CartesianVector &point : trajPoints)
        {
            mcParticleTrajX.push_back(point.GetX());
            mcParticleTrajY.push_back(point.GetY());
            mcParticleTrajZ.push_back(point.GetZ());
            mcParticleTrajT.push_back(pLArMC->GetT0());            
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitPEs", &opHitPEs));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitPeakTimes", &opHitPeakTimes));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitStartTimes", &opHitStartTimes));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitWidths", &opHitWidths));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitPositionsX", &opHitPositionsX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitPositionsY", &opHitPositionsY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitPositionsZ", &opHitPositionsZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "opHitChannels", &opHitChannels));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpTimes", &mcParticleTimes));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpEnergies", &mcParticleEnergies));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpPdgs", &mcParticlePdgs));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpIsMichel", &mcParticleIsMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpTrajPointX", &mcParticleTrajX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpTrajPointY", &mcParticleTrajY));        
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpTrajPointZ", &mcParticleTrajZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mcpTrajPointT", &mcParticleTrajT));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SolutionLightAlgorithm::IsMichelElectron(const LArMCParticle *const pMC) const
{
    if (pMC->GetParticleId() != E_MINUS)
        return false;

    if (pMC->GetParentList().size() != 1 || pMC->GetParentList().front()->GetParticleId() != MU_MINUS)
        return false;

    // This should be the last step but the G4 process flag is bugged to always be capture for mu-...
    // if (pLArMC->GetProcess() != MC_PROC_DECAY)
    //     return false

    // Sub MeV and prompt -> probably the electron ejected when the muon comes to rest in an orbital of an Ar atom... really annoying the process info is broken
    if (pMC->GetT0() < 10 && pMC->GetEnergy() < 0.001)
        return false;

    const MCParticle *const pParentMuonMC{pMC->GetParentList().front()};
    auto sameVertex = [](const MCParticle *const pMC1, const MCParticle *const pMC2, const float eps=0.001) -> bool
      { return std::abs(pMC1->GetVertex().GetX() - pMC2->GetVertex().GetX()) < eps &&
               std::abs(pMC1->GetVertex().GetY() - pMC2->GetVertex().GetY()) < eps &&
               std::abs(pMC1->GetVertex().GetZ() - pMC2->GetVertex().GetZ()) < eps; };
    bool foundSiblingNumu{false}, foundSiblingNuebar{false};
    for (const MCParticle *const pSiblingMC : pParentMuonMC->GetDaughterList())
    {
        if (pSiblingMC == pMC)
            continue;

        if (!sameVertex(pSiblingMC, pMC))
            continue;

        if (pSiblingMC->GetParticleId() == NU_MU)
            foundSiblingNumu = true;
        else if (pSiblingMC->GetParticleId() == NU_E_BAR)
            foundSiblingNuebar = true;
    }
    if (!foundSiblingNumu || !foundSiblingNuebar)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionLightAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OpticalHitListName", m_opticalHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

