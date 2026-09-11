/**
 *  @file   larpandoracontent/LArSolutions/SolutionLightClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the SolutionLightClusteringAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArSolutions/SolutionLightClusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SolutionLightClusteringAlgorithm::SolutionLightClusteringAlgorithm() :
    m_opHitPEThreshold(1000.f),
    m_maxDeltaY(100.f),
    m_maxDeltaZ(100.f),
    m_maxDeltaT(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionLightClusteringAlgorithm::Run()
{
    const CaloHitList *pOpHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_opticalHitListName, pOpHitList));

    // Filter and sort list for reproducibility
    OpHitVector opHitVector;
    this->PrepareOpHitList(pOpHitList, opHitVector);
    
    // Create temporary list to store new clusters
    const ClusterList *pTemporaryList(nullptr);
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTemporaryList, temporaryListName));

    // Make clusters   
    bool seedFound(true);
    OpHitVector consideredSeeds;
    
    while(seedFound)
    {
        const LArOpHit *pOpSeed(nullptr);          
        seedFound = this->GetSeed(opHitVector, consideredSeeds, pOpSeed);

        if (!seedFound)
            break;
        
        consideredSeeds.push_back(pOpSeed);
        
        CaloHitList opClusterHits;
        this->GetClusterHits(pOpSeed, opHitVector, opClusterHits);
        
        // Create cluster with collected hits
        const Cluster *pCluster(nullptr);
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList.insert(parameters.m_caloHitList.begin(), opClusterHits.begin(), opClusterHits.end());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
    }

    if (!pTemporaryList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputOpClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputOpClusterListName));
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionLightClusteringAlgorithm::PrepareOpHitList(const CaloHitList *const pOpHitList, OpHitVector &opHitVector) const
{
    for (const CaloHit *const pCaloHit : *pOpHitList)
    {
        const LArOpHit *const pOpHit(dynamic_cast<const LArOpHit *>(pCaloHit));
        
        if (!pOpHit)
            continue;

        if (pOpHit->GetInputEnergy() < m_opHitPEThreshold)
            continue;

        opHitVector.push_back(pOpHit);
    }
        
    // Sort list for reproducibility
    std::sort(opHitVector.begin(), opHitVector.end(), [](const LArOpHit *const pHit1, const LArOpHit *const pHit2)
    {
        return pHit1->GetInputEnergy() > pHit2->GetInputEnergy();
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SolutionLightClusteringAlgorithm::GetSeed(const OpHitVector &opHitVector, const OpHitVector &consideredSeeds, const LArOpHit *& pOpSeed) const
{
    bool seedFound(false);
    float highestPE(0.f);
    
    for (const LArOpHit *const pOpHit : opHitVector)
    {
        if (std::find(consideredSeeds.begin(), consideredSeeds.end(), pOpHit) != consideredSeeds.end())
            continue;
        
        const CaloHit *const pCaloHit(pOpHit);
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        if (pOpHit->GetInputEnergy() > highestPE)
        {
            highestPE = pOpHit->GetInputEnergy();
            pOpSeed = pOpHit;
            seedFound = true;
        }
    }

    return seedFound;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionLightClusteringAlgorithm::GetClusterHits(const LArOpHit *const pOpSeed, OpHitVector &opHitVector, CaloHitList &opClusterHits) const
{
    opClusterHits.push_back(pOpSeed);

    bool clusterModified(true);

    while(clusterModified)
    {
        clusterModified = false;
        
        for (const LArOpHit *const pEventOpHit : opHitVector)
        {            
            const CaloHit *const pEventCaloHit(pEventOpHit);
            if (!PandoraContentApi::IsAvailable(*this, pEventCaloHit))
                continue;

            // Make sure that we haven't already collected this hit
            if (std::find(opClusterHits.begin(), opClusterHits.end(), pEventCaloHit) != opClusterHits.end())
                continue;

            bool isClose(false);
            for (const CaloHit *const pClusterCaloHit : opClusterHits)
            {
                const LArOpHit *const pClusterOpHit(dynamic_cast<const LArOpHit *>(pClusterCaloHit));
       
                if (!pClusterOpHit)
                    continue;
                
                if (!this->IsClose(pEventOpHit, pClusterOpHit))
                    continue;

                isClose = true;
                clusterModified = true;
                break;
            }

            if (isClose)
                opClusterHits.push_back(pEventOpHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SolutionLightClusteringAlgorithm::IsClose(const LArOpHit *const pEventOpHit, const LArOpHit *const pClusterOpHit) const
{
    const float deltaY(std::fabs(pEventOpHit->GetPositionVector().GetY() - pClusterOpHit->GetPositionVector().GetY()));
    const float deltaZ(std::fabs(pEventOpHit->GetPositionVector().GetZ() - pClusterOpHit->GetPositionVector().GetZ()));
    const float deltaT(std::fabs(pEventOpHit->GetTime() - pClusterOpHit->GetTime()));
            
    return ((deltaY < m_maxDeltaY) && (deltaZ < m_maxDeltaZ) && (deltaT < m_maxDeltaT));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionLightClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OpticalHitListName", m_opticalHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputOpClusterListName", m_outputOpClusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "OpHitPEThreshold", m_opHitPEThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDeltaY", m_maxDeltaY));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDeltaZ", m_maxDeltaZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDeltaT", m_maxDeltaT));        
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

