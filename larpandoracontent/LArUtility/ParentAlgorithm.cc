/**
 *  @file   larpandoracontent/LArUtility/ParentAlgorithm.cc
 *
 *  @brief  Implementation of the parent algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArUtility/ParentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ParentAlgorithm::ParentAlgorithm() :
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldIdentifyNeutrinoSlice(true),
    m_printStatus(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentAlgorithm::Run()
{
    if (m_shouldRunAllHitsCosmicReco)
    {
        PfoList parentCosmicRayPfos;
        this->RunAllHitsCosmicRayReconstruction(parentCosmicRayPfos);

        if (parentCosmicRayPfos.empty())
            return STATUS_CODE_SUCCESS;

        if (m_shouldRunCosmicHitRemoval)
        {
            PfoList ambiguousPfos;
            m_pCosmicRayTaggingTool->FindAmbiguousPfos(parentCosmicRayPfos, ambiguousPfos);
            this->RemoveAmbiguousCosmicRayPfos(ambiguousPfos);
        }

        this->RunAlgorithm(m_crListMovingAlgorithm);
    }

    if (!m_shouldRunNeutrinoRecoOption && !m_shouldRunCosmicRecoOption)
        return STATUS_CODE_SUCCESS;

    SliceList sliceList;
    this->RunFastReconstructionAndSlicing(sliceList);

    if (sliceList.empty())
        return STATUS_CODE_SUCCESS;

    SliceIndexToPfoListMap sliceToCosmicRayPfosMap;
    SliceIndexToPropertiesMap sliceIndexToPropertiesMap;
    this->ReconstructSlices(sliceList, sliceToCosmicRayPfosMap, sliceIndexToPropertiesMap);

    if (sliceToCosmicRayPfosMap.empty())
        return STATUS_CODE_SUCCESS;

    unsigned int neutrinoSliceIndex(0);
    if (m_shouldIdentifyNeutrinoSlice && m_pNeutrinoIdTool->GetNeutrinoSliceIndex(sliceIndexToPropertiesMap, neutrinoSliceIndex))
    {
        this->RemoveSliceCosmicRayReconstruction(sliceToCosmicRayPfosMap, neutrinoSliceIndex);
        this->AddSliceNeutrinoReconstruction(sliceList, neutrinoSliceIndex);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RunAllHitsCosmicRayReconstruction(PfoList &parentCosmicRayPfos) const
{
    SliceList sliceList;
    this->CopyAllHitsToSingleSlice(sliceList);
    this->CosmicRayReconstruction(sliceList.front(), std::string());

    const PfoList *pParentCRPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_crParentListName, pParentCRPfoList));

    if (pParentCRPfoList)
        parentCosmicRayPfos = *pParentCRPfoList;

    if (m_printStatus || PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: cosmic-ray reconstuction done" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RemoveAmbiguousCosmicRayPfos(const PfoList &ambiguousPfos) const
{
    const PfoList *pParentCosmicRayPfos(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_crParentListName, pParentCosmicRayPfos));

    const PfoList *pCosmicRayDaughterPfos(nullptr);
    (void) PandoraContentApi::GetList(*this, m_crDaughterListName, pCosmicRayDaughterPfos);

    for (const Pfo *const pPfo : ambiguousPfos)
    {
        const std::string listName(std::find(pParentCosmicRayPfos->begin(), pParentCosmicRayPfos->end(), pPfo) != pParentCosmicRayPfos->end() ? m_crParentListName :
            (pCosmicRayDaughterPfos && std::find(pCosmicRayDaughterPfos->begin(), pCosmicRayDaughterPfos->end(), pPfo) != pCosmicRayDaughterPfos->end()) ? m_crDaughterListName : "");

        if (listName.empty())
        {
            std::cout << "ParentAlgorithm::RemoveAmbiguousCosmicRayPfos: parent/daughter particles list mismatch" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, listName));
    }

    this->RunAlgorithm(m_crListPruningAlgorithm);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RunFastReconstructionAndSlicing(SliceList &sliceList) const
{
    if (m_shouldRunSlicing)
    {
        this->RunSlicing(sliceList);
    }
    else
    {
        this->CopyAllHitsToSingleSlice(sliceList);
    }

    if (m_printStatus || PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: slicing done, nSlices " << sliceList.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::ReconstructSlices(const SliceList &sliceList, SliceIndexToPfoListMap &sliceToCosmicRayPfosMap,
    SliceIndexToPropertiesMap &sliceIndexToPropertiesMap) const
{
    unsigned int sliceIndex(0);

    for (const Slice &slice : sliceList)
    {
        const unsigned int thisSliceIndex(sliceIndex++);
        const std::string sliceIndexString(TypeToString(thisSliceIndex));
        SliceProperties sliceProperties;

        if (m_shouldRunNeutrinoRecoOption)
        {
            this->NeutrinoReconstruction(slice, sliceIndexString);

            const PfoList *pParentNeutrinosInSlice(nullptr);
            (void) PandoraContentApi::GetList(*this, m_nuParentListName, pParentNeutrinosInSlice);

            if (pParentNeutrinosInSlice)
                m_pNeutrinoIdTool->FillNeutrinoProperties(pParentNeutrinosInSlice, sliceProperties);

            const std::string postProcessAlg(!m_shouldRunCosmicRecoOption ? m_nuListMovingAlgorithm : m_nuListDeletionAlgorithm);
            this->RunAlgorithm(postProcessAlg);

            if (m_printStatus || PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ParentAlgorithm: neutrino reconstruction done for slice " << thisSliceIndex << std::endl;
        }

        if (m_shouldRunCosmicRecoOption)
        {
            this->CosmicRayReconstruction(slice, sliceIndexString);

            const PfoList *pParentCRPfosInSlice(nullptr);
            (void) PandoraContentApi::GetList(*this, m_crParentListName, pParentCRPfosInSlice);

            if (pParentCRPfosInSlice)
                m_pNeutrinoIdTool->FillCosmicRayProperties(pParentCRPfosInSlice, sliceProperties);

            if (pParentCRPfosInSlice)
                sliceToCosmicRayPfosMap[thisSliceIndex] = *pParentCRPfosInSlice;

            this->RunAlgorithm(m_crListMovingAlgorithm);

            if (m_printStatus || PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ParentAlgorithm: cosmic-ray reconstruction done for slice " << thisSliceIndex << std::endl;
        }

        sliceIndexToPropertiesMap[thisSliceIndex] = sliceProperties;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::RemoveSliceCosmicRayReconstruction(const SliceIndexToPfoListMap &sliceToCosmicRayPfosMap, const unsigned int neutrinoSliceIndex) const
{
    const PfoList *pParentCosmicRayPfos(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_outputListPrefix + m_crParentListName, pParentCosmicRayPfos));

    const PfoList *pCosmicRayDaughterPfos(nullptr);
    (void) PandoraContentApi::GetList(*this, m_outputListPrefix + m_crDaughterListName, pCosmicRayDaughterPfos);

    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(sliceToCosmicRayPfosMap.at(neutrinoSliceIndex), allPfosToDelete);

    for (const Pfo *const pPfo : allPfosToDelete)
    {
        const std::string listName(std::find(pParentCosmicRayPfos->begin(), pParentCosmicRayPfos->end(), pPfo) != pParentCosmicRayPfos->end() ? m_outputListPrefix + m_crParentListName :
            (pCosmicRayDaughterPfos && std::find(pCosmicRayDaughterPfos->begin(), pCosmicRayDaughterPfos->end(), pPfo) != pCosmicRayDaughterPfos->end()) ? m_outputListPrefix + m_crDaughterListName : "");

        if (listName.empty())
        {
            std::cout << "ParentAlgorithm::RemoveSliceCosmicRayReconstruction: parent/daughter particles list mismatch" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, listName));
    }

    this->RunAlgorithm(m_outputListPruningAlgorithm);

    if (m_printStatus || PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: cosmic-ray reconstruction removed for slice " << neutrinoSliceIndex << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::AddSliceNeutrinoReconstruction(const SliceList &sliceList, const unsigned int neutrinoSliceIndex) const
{
    this->NeutrinoReconstruction(sliceList.at(neutrinoSliceIndex), TypeToString(neutrinoSliceIndex));
    this->RunAlgorithm(m_nuListMovingAlgorithm);

    if (m_printStatus || PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "ParentAlgorithm: neutrino reconstruction applied to slice " << neutrinoSliceIndex << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::FastReconstruction() const
{
    this->RunTwoDClustering(std::string(), m_nuClusteringAlgorithm, false, m_nuTwoDAlgorithms);
    this->RunAlgorithms(m_nuThreeDAlgorithms);
    this->RunAlgorithms(m_nuThreeDHitAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::CosmicRayReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    this->SaveTwoDCaloHitLists(slice, sliceIndexString);
    this->RunTwoDClustering(sliceIndexString, m_crTrackClusteringAlgorithm, false, m_crTwoDAlgorithms);
    this->RunAlgorithms(m_crThreeDAlgorithms);
    this->RunAlgorithm(m_crListPruningAlgorithm);
    this->RunTwoDClustering(sliceIndexString, m_crDeltaRayClusteringAlgorithm, true);
    this->RunAlgorithms(m_crDeltaRayAlgorithms);

    for (const HitType hitType : m_hitTypeList)
    {
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, (PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType))));
        this->RunAlgorithms(m_crTwoDRemnantAlgorithms);
    }

    this->RunAlgorithms(m_crThreeDRemnantAlgorithms);
    this->RunAlgorithms(m_crThreeDHitAlgorithms);
    this->RunAlgorithms(m_crVertexAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentAlgorithm::NeutrinoReconstruction(const ParentSlicingBaseAlgorithm::Slice &slice, const std::string &sliceIndexString) const
{
    this->SaveTwoDCaloHitLists(slice, sliceIndexString);
    this->RunTwoDClustering(sliceIndexString, m_nuClusteringAlgorithm, false, m_nuTwoDAlgorithms);
    this->RunAlgorithms(m_nuVertexAlgorithms);
    this->RunAlgorithms(m_nuThreeDAlgorithms);
    this->RunAlgorithms(m_nuTwoDMopUpAlgorithms);
    this->RunAlgorithms(m_nuThreeDHitAlgorithms);
    this->RunAlgorithms(m_nuThreeDMopUpAlgorithms);
    this->RunAlgorithms(m_nuNeutrinoAlgorithms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldRunAllHitsCosmicReco", m_shouldRunAllHitsCosmicReco));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldRunCosmicHitRemoval", m_shouldRunCosmicHitRemoval));

    if (m_shouldRunCosmicHitRemoval && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "ParentAlgorithm::ReadSettings - ShouldRunCosmicHitRemoval requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
            "CosmicRayTagging", pAlgorithmTool));

        m_pCosmicRayTaggingTool = dynamic_cast<CosmicRayTaggingBaseTool*>(pAlgorithmTool);

        if (!m_pCosmicRayTaggingTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldRunNeutrinoRecoOption", m_shouldRunNeutrinoRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldRunCosmicRecoOption", m_shouldRunCosmicRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldIdentifyNeutrinoSlice", m_shouldIdentifyNeutrinoSlice));

    if (m_shouldIdentifyNeutrinoSlice && (!m_shouldRunSlicing || !m_shouldRunNeutrinoRecoOption || !m_shouldRunCosmicRecoOption))
    {
        std::cout << "ParentAlgorithm::ReadSettings - ShouldIdentifyNeutrinoSlice requires ShouldRunSlicing and both neutrino and cosmic reconstruction options" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldIdentifyNeutrinoSlice)
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
            "NeutrinoId", pAlgorithmTool));

        m_pNeutrinoIdTool = dynamic_cast<NeutrinoIdBaseTool*>(pAlgorithmTool);

        if (!m_pNeutrinoIdTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintStatus", m_printStatus));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CRParentListName", m_crParentListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CRDaughterListName", m_crDaughterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NuParentListName", m_nuParentListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputListPrefix", m_outputListPrefix));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "OutputListPruning", m_outputListPruningAlgorithm));

    TiXmlElement *const pCosmicRayXmlElement(xmlHandle.FirstChild("CosmicRayReconstruction").Element());

    if ((m_shouldRunAllHitsCosmicReco || m_shouldRunCosmicRecoOption) && pCosmicRayXmlElement)
    {
        const TiXmlHandle crXmlHandle(pCosmicRayXmlElement);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "TwoDTrackClustering", m_crTrackClusteringAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "TwoDDeltaRayClustering", m_crDeltaRayClusteringAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "TwoDAlgorithms", m_crTwoDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "ThreeDAlgorithms", m_crThreeDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "DeltaRayAlgorithms", m_crDeltaRayAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "TwoDRemnantAlgorithms", m_crTwoDRemnantAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "ThreeDRemnantAlgorithms", m_crThreeDRemnantAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "ThreeDHitAlgorithms", m_crThreeDHitAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, crXmlHandle,
            "VertexAlgorithms", m_crVertexAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "ListPruning", m_crListPruningAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "ListDeletion", m_crListDeletionAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, crXmlHandle,
            "ListMoving", m_crListMovingAlgorithm));
    }

    TiXmlElement *const pNeutrinoXmlElement(xmlHandle.FirstChild("NeutrinoReconstruction").Element());

    if ((m_shouldRunSlicing || m_shouldRunNeutrinoRecoOption) && pNeutrinoXmlElement)
    {
        const TiXmlHandle nuXmlHandle(pNeutrinoXmlElement);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, nuXmlHandle,
            "TwoDClustering", m_nuClusteringAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "TwoDAlgorithms", m_nuTwoDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "ThreeDAlgorithms", m_nuThreeDAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "ThreeDHitAlgorithms", m_nuThreeDHitAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "VertexAlgorithms", m_nuVertexAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "TwoDMopUpAlgorithms", m_nuTwoDMopUpAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "ThreeDMopUpAlgorithms", m_nuThreeDMopUpAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, nuXmlHandle,
            "NeutrinoAlgorithms", m_nuNeutrinoAlgorithms));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, nuXmlHandle,
            "ListDeletion", m_nuListDeletionAlgorithm));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, nuXmlHandle,
            "ListMoving", m_nuListMovingAlgorithm));
    }

    return ParentSlicingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
