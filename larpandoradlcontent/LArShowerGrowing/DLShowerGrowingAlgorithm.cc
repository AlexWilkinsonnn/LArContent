/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning shower growing algorithm
 *
 *  $Log: $
 */

#include "larpandoradlcontent/LArShowerGrowing/DLShowerGrowingAlgorithm.h"

using namespace pandora;
// using namespace lar_content;

namespace lar_dl_content
{

DLShowerGrowingAlgorithm::DLShowerGrowingAlgorithm() :
    m_trainingMode{false},
    m_trainingTreeName{""}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLShowerGrowingAlgorithm::~DLShowerGrowingAlgorithm()
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLShowerGrowingAlgorithm::Run()
{
    if (m_trainingMode)
    {
        return this->PrepareTrainingSample();
    }
    else
    {
        return this->Infer();
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLShowerGrowingAlgorithm::PrepareTrainingSample()
{
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLShowerGrowingAlgorithm::Infer()
{
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingTreeName", m_trainingTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
