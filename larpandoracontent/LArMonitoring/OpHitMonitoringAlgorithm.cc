/**
 *  @file   larpandoracontent/LArMonitoring/OpHitMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the optical hit monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/OpHitMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

using namespace pandora;

namespace lar_content
{

OpHitMonitoringAlgorithm::OpHitMonitoringAlgorithm() :
    m_opticalTimeExtent{50.f},
    m_opticalMagnitudeScale{0.1f},
    m_opticalTimeMin{0.f},
    m_opticalTimeMax{20.f},
    m_simpleMode{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

OpHitMonitoringAlgorithm::~OpHitMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OpHitMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    for (const auto pfoListName : m_pfoListNames)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->VisualizePfos(pfoListName));
    }

    const CaloHitList *pOpHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_opticalHitListName, pOpHitList));

    std::map<unsigned int, CartesianVector> detectorPositions;

    for (const CaloHit *const pCaloHit : *pOpHitList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->VisualizeOpHit(pCaloHit, detectorPositions));
    }

    if (!m_simpleMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->VisualizeOpDetTimeRanges(detectorPositions));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OpHitMonitoringAlgorithm::VisualizePfos(const std::string &listName) const
{
    const PfoList *pPfoList{nullptr};
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pPfoList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "OpHitMonitoringAlgorithm: pfo list " << listName << " unavailable." << std::endl;
        return STATUS_CODE_SUCCESS;
    }
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, listName.c_str(), AUTOITER, true, true));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OpHitMonitoringAlgorithm::VisualizeOpHit(const CaloHit *const pCaloHit,
    std::map<unsigned int, CartesianVector> &detectorPositions) const
{
    const LArOpHit *const pOpHit{dynamic_cast<const LArOpHit *>(pCaloHit)};
    if (!pOpHit)
        return STATUS_CODE_FAILURE;

    const CartesianVector pos{pOpHit->GetPositionVector()};
    const float magnitudeLength{pOpHit->GetInputEnergy() * m_opticalMagnitudeScale};

    if (m_simpleMode)
    {
        const int markerSize{static_cast<int>(std::ceil(magnitudeLength))};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "Optical Hit", ORANGE, markerSize));
        return STATUS_CODE_SUCCESS;
    }

    const float peakTime{pOpHit->GetTime()};
    const float displayTime{std::max(m_opticalTimeMin, std::min(peakTime, m_opticalTimeMax))};

    const float halfWidth{0.5f * pOpHit->GetWidth()};
    const float xStart{pos.GetX() - (0.5f * m_opticalTimeExtent) + (displayTime - m_opticalTimeMin - halfWidth) * m_opticalTimeScale};
    const float xEnd{pos.GetX() - (0.5f * m_opticalTimeExtent) + (displayTime - m_opticalTimeMin + halfWidth) * m_opticalTimeScale};
    const CartesianVector timeXStartPos(xStart, pos.GetY(), pos.GetZ());
    const CartesianVector timeXEndPos(xEnd, pos.GetY(), pos.GetZ());

    const CartesianVector peakPos((xStart + xEnd) * 0.5f, pos.GetY(), pos.GetZ());
    const CartesianVector magnitudeYStartPos(peakPos.GetX(), peakPos.GetY() - magnitudeLength, peakPos.GetZ());
    const CartesianVector magnitudeYEndPos(peakPos.GetX(), peakPos.GetY() + magnitudeLength, peakPos.GetZ());
    const CartesianVector magnitudeZStartPos(peakPos.GetX(), peakPos.GetY(), peakPos.GetZ() - magnitudeLength);
    const CartesianVector magnitudeZEndPos(peakPos.GetX(), peakPos.GetY(), peakPos.GetZ() + magnitudeLength);

    const Color hitColor{(peakTime < m_opticalTimeMin) || (peakTime > m_opticalTimeMax) ? RED : ORANGE}; // If overflow/underflow time

    detectorPositions.emplace(pOpHit->GetChannel(), pos);
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &timeXStartPos, &timeXEndPos, "OpHit time width", hitColor, 6, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(
        this->GetPandora(), &magnitudeYStartPos, &magnitudeYEndPos, "OpHit magnitude Y", hitColor, 6, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(
        this->GetPandora(), &magnitudeZStartPos, &magnitudeZEndPos, "OpHit magnitude Z", hitColor, 6, 1));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OpHitMonitoringAlgorithm::VisualizeOpDetTimeRanges(const std::map<unsigned int, CartesianVector> &detectorPositions) const
{
    for (const auto &detectorEntry : detectorPositions)
    {
        const CartesianVector &pos{detectorEntry.second};
        const CartesianVector startPos(pos.GetX() - (0.5f * m_opticalTimeExtent), pos.GetY(), pos.GetZ());
        const CartesianVector endPos(pos.GetX() + (0.5f * m_opticalTimeExtent), pos.GetY(), pos.GetZ());
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPos, &endPos, "Optical detector time range", BLUE, 1, 1));
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OpHitMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OpticalHitListName", m_opticalHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF( STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OpticalTimeExtent", m_opticalTimeExtent));
    if (m_opticalTimeExtent <= 0.f)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OpticalMagnitudeScale", m_opticalMagnitudeScale));
    if (m_opticalMagnitudeScale <= 0.f)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OpticalTimeMin", m_opticalTimeMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OpticalTimeMax", m_opticalTimeMax));
    if (m_opticalTimeMax <= m_opticalTimeMin)
        return STATUS_CODE_INVALID_PARAMETER;
    m_opticalTimeScale = m_opticalTimeExtent / (m_opticalTimeMax - m_opticalTimeMin);

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SimpleMode", m_simpleMode));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
