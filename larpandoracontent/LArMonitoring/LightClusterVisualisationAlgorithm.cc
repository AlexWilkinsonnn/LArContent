/**
 *  @file   larpandoracontent/LArMonitoring/LightClusterVisualisationAlgorithm.cc
 *
 *  @brief  Implementation of the LightClusterVisualisationAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArMonitoring/LightClusterVisualisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

LightClusterVisualisationAlgorithm::LightClusterVisualisationAlgorithm() :
    m_opticalMagnitudeScale(0.001f),
    m_minT(0.f),
    m_maxT(10.f),
    m_minClusterHits(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LightClusterVisualisationAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));

    const ClusterList *pOpClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_opClusterListName, pOpClusterList));
    ClusterVector pOpClusters(pOpClusterList->begin(), pOpClusterList->end());

    // Sort wrt time
    std::sort(pOpClusters.begin(), pOpClusters.end(), [this](const Cluster *const pCluster1, const Cluster *const pCluster2)
    {
        const float meanTime1(this->GetClusterT0(pCluster1));
        const float meanTime2(this->GetClusterT0(pCluster2));            
        
        return meanTime1 < meanTime2;
    });

    for (const Cluster *const pOpCluster : pOpClusters)
    {
        const float clusterT0(this->GetClusterT0(pOpCluster));

        if ((clusterT0 < m_minT) || (clusterT0 > m_maxT))
            continue;

        CaloHitList opCaloHits;
        LArClusterHelper::GetAllHits(pOpCluster, opCaloHits);

        if (opCaloHits.size() < m_minClusterHits)
            continue;
        
        std::cout << "---------------------------" << std::endl;
        std::cout << "Cluster T0: " << clusterT0  << std::endl;
        std::cout << "---------------------------" << std::endl;

        for (const CaloHit *const pOpCaloHit : opCaloHits)
        {
            const LArOpHit *const pOpHit(dynamic_cast<const LArOpHit *>(pOpCaloHit));
        
            if (!pOpHit)
                return STATUS_CODE_FAILURE;

            const CartesianVector pos(pOpHit->GetPositionVector());
            const float magnitudeLength(pOpHit->GetInputEnergy() * m_opticalMagnitudeScale);
            const int markerSize(static_cast<int>(std::ceil(magnitudeLength)));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, ("Optical Hit - " + std::to_string(pOpHit->GetInputEnergy())), ORANGE, markerSize));
        }
        
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
    
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "Input", pMCParticleList));

    for (const MCParticle *const pMC : *pMCParticleList)
    {
        const LArMCParticle *const pLArMC(dynamic_cast<const LArMCParticle *>(pMC));
        
        if (!pLArMC)
            return STATUS_CODE_INVALID_PARAMETER;
        
        CartesianPointVector trajPoints(pLArMC->GetTrajPoints());
        for (const CartesianVector &point : trajPoints)
        {
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &point, "MC Trajectory", BLUE, 2));
        }        
        
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LightClusterVisualisationAlgorithm::GetClusterT0(const Cluster *const pOpCluster)
{
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pOpCluster, clusterHits);

    float timeSum(0.f);

    for (const CaloHit *const pCaloHit : clusterHits)
    {
        const LArOpHit *const pOpHit(dynamic_cast<const LArOpHit *>(pCaloHit));
        
        if (!pOpHit)
            continue;
        timeSum += pOpHit->GetTime();
    }

    return (timeSum / clusterHits.size());
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LightClusterVisualisationAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OpClusterListName", m_opClusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OpticalMagnitudeScale", m_opticalMagnitudeScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinT", m_minT));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinT", m_minT));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

