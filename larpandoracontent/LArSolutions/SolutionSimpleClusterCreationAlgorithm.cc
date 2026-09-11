/**
 *  @file   larpandoracontent/LArSolutions/SolutionSimpleClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the SolutionSimpleClusterCreationAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArSolutions/SolutionSimpleClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SolutionSimpleClusterCreationAlgorithm::SolutionSimpleClusterCreationAlgorithm() :
    m_maxNCaloHits(5),
    m_maxSeparation(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionSimpleClusterCreationAlgorithm::Run()
{
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // Sort list for reproducibility
    CaloHitVector caloHitVector(pCaloHitList->begin(), pCaloHitList->end()); // Remember that our pandora list is const!
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    // Create temporary list to store new clusters
    const ClusterList *pTemporaryList(nullptr);
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTemporaryList, temporaryListName));

    // Now form clusters
    const Cluster *pCluster(nullptr);
    
    for (const CaloHit *const pSeedHit : caloHitVector)
    {
        if (!PandoraContentApi::IsAvailable(*this, pSeedHit))
            continue;

        CaloHitList newClusterHits;
        this->GetClusterHits(pSeedHit, caloHitVector, newClusterHits);
        
        // Create cluster with collected hits
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList.insert(parameters.m_caloHitList.begin(), newClusterHits.begin(), newClusterHits.end());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

        // Visualise
        //CartesianVector seedPosition(pSeedHit->GetPositionVector());
        //ClusterList visClusters({pCluster});
        //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &seedPosition, "Seed", VIOLET, 1.5));
        //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visClusters, "New Cluster", BLACK));
    }

    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    if (!pTemporaryList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputClusterListName));
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionSimpleClusterCreationAlgorithm::GetClusterHits(const CaloHit *const pSeedHit, const CaloHitVector &eventHits, CaloHitList &clusterHits) const
{
    // Add seed hit to hit collection
    clusterHits.push_back(pSeedHit);

    // Grow collection of hits until we can't find a hit close enough to add, or the group becomes too large        
    bool found(true);
    
    while (found && (clusterHits.size() < m_maxNCaloHits))
    {
        found = false;

        // Find the closest event hit to add to the group
        float smallestSep(std::numeric_limits<float>::max());
        const CaloHit *pBestHit(nullptr);
   
        for (const CaloHit *const pEventHit : eventHits)
        {
            if (!PandoraContentApi::IsAvailable(*this, pEventHit))
                continue;

            // Make sure that we haven't already collected this hit
            if (std::find(clusterHits.begin(), clusterHits.end(), pEventHit) != clusterHits.end())
                continue;
            
            const float separation(LArClusterHelper::GetClosestDistance(pEventHit->GetPositionVector(), clusterHits));

            // Move on if the hit is too far away
            if (separation > m_maxSeparation)
                continue;

            if (separation < smallestSep)
            {
                found = true; 
                smallestSep = separation;
                pBestHit = pEventHit;
            }
        }

        if (found && pBestHit)
            clusterHits.push_back(pBestHit);            
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionSimpleClusterCreationAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxNCaloHits", m_maxNCaloHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeparation", m_maxSeparation));
    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

