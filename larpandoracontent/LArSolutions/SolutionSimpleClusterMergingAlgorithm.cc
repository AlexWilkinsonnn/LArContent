/**
 *  @file   larpandoracontent/LArSolutions/SolutionSimpleClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the SolutionSimpleClusterMergingAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArSolutions/SolutionSimpleClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SolutionSimpleClusterMergingAlgorithm::SolutionSimpleClusterMergingAlgorithm() :
    m_minClusterHits(2),
    m_maxClusterSeparation(3.f),
    m_maxOpeningAngle(25.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionSimpleClusterMergingAlgorithm::Run()
{
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName));

    ClusterVector filteredClusters;
    this->GetFilteredClusters(pClusterList, filteredClusters);

    // Merge clusters!
    std::set<const Cluster*> staleClusters;
    for (const Cluster *const pParentCluster : filteredClusters)
    {        
        if (staleClusters.count(pParentCluster))
            continue;

        for (const Cluster *const pChildCluster : filteredClusters)
        {
            if ((pParentCluster == pChildCluster) || staleClusters.count(pChildCluster))
                continue;
            
            if (!this->AreClustersAssociated(pParentCluster, pChildCluster))
                continue;
            
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pChildCluster));
            staleClusters.insert(pChildCluster);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionSimpleClusterMergingAlgorithm::GetFilteredClusters(const ClusterList *const pClusterList, ClusterVector &filteredClusters) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minClusterHits)
            continue;
        
        filteredClusters.push_back(pCluster);
    }

    std::sort(filteredClusters.begin(), filteredClusters.end(), LArClusterHelper::SortByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SolutionSimpleClusterMergingAlgorithm::AreClustersAssociated(const Cluster *const pParentCluster, const Cluster *const pChildCluster) const
{
    CartesianVector parentInnerCoord(0.f, 0.f, 0.f);
    CartesianVector parentOuterCoord(0.f, 0.f, 0.f);
    CartesianVector childInnerCoord(0.f, 0.f, 0.f);
    CartesianVector childOuterCoord(0.f, 0.f, 0.f);
    
    LArClusterHelper::GetExtremalCoordinates(pParentCluster, parentInnerCoord, parentOuterCoord);
    LArClusterHelper::GetExtremalCoordinates(pChildCluster, childInnerCoord, childOuterCoord);

    float separation(std::numeric_limits<float>::max());
    for (const CartesianVector &parentCoord : {parentInnerCoord, parentOuterCoord})
    {
        for (const CartesianVector &childCoord : {childInnerCoord, childOuterCoord})
        {
            const float thisSeparation((parentCoord - childCoord).GetMagnitude());
            separation = std::min(separation, thisSeparation);
        }
    }

    if (separation > m_maxClusterSeparation)
        return false;

    const CartesianVector parentDirection((parentInnerCoord - parentOuterCoord).GetUnitVector());
    const CartesianVector childDirection((childInnerCoord - childOuterCoord).GetUnitVector());
    const float openingAngle(parentDirection.GetOpeningAngle(childDirection) * 180.f / M_PI);

    if (openingAngle > m_maxOpeningAngle)
        return false;

    // Monitoring
    // ClusterList visParentCluster({pParentCluster});
    // ClusterList visChildCluster({pChildCluster});
    // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visParentCluster, "Parent", BLACK));
    // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visChildCluster, "Child", VIOLET));
    // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &parentInnerCoord, "Parent Inner", SPRING, 1.5));
    // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &parentOuterCoord, "Parent Outer", SPRING, 1.5));
    // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &childInnerCoord, "Child Inner", SPRING, 1.5));
    // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &childOuterCoord, "Child Outer", SPRING, 1.5));
    // const CartesianVector parentDir1(((parentInnerCoord + parentOuterCoord) * 0.5f) + (parentDirection * 5.f));
    // const CartesianVector parentDir2(((parentInnerCoord + parentOuterCoord) * 0.5f) - (parentDirection * 5.f));
    // PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &parentDir1, &parentDir2, "Parent Direction", BLACK, 2, 1));
    // const CartesianVector childDir1(((childInnerCoord + childOuterCoord) * 0.5f) + (childDirection * 5.f));
    // const CartesianVector childDir2(((childInnerCoord + childOuterCoord) * 0.5f) - (childDirection * 5.f));    
    // PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &childDir1, &childDir2, "Child Direction", VIOLET, 2, 1));
    // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    
    return true;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionSimpleClusterMergingAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterSeparation", m_maxClusterSeparation));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxOpeningAngle", m_maxOpeningAngle));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

