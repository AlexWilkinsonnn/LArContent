/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/SolutionSimpleClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the SolutionSimpleClusterSplittingAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArSolutions/SolutionSimpleClusterSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SolutionSimpleClusterSplittingAlgorithm::SolutionSimpleClusterSplittingAlgorithm() :
    m_nAngularBins(20),
    m_radiusForDirEstimate(5.f),
    m_maxDistToSplitPos(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionSimpleClusterSplittingAlgorithm::Run()
{
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));

    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName));
    
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const MCParticle *pMCMuon(nullptr);
    this->GetMCMuon(pMCParticleList, pMCMuon);
    
    const CartesianVector muonVertex(pMCMuon->GetVertex());
    const CartesianVector vertexProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), muonVertex, LArClusterHelper::GetClusterHitType(pClusterList->front())));
    const CartesianVector muonEndpoint(pMCMuon->GetEndpoint());
    const CartesianVector endpointProjection(LArGeometryHelper::ProjectPosition(this->GetPandora(), muonEndpoint, LArClusterHelper::GetClusterHitType(pClusterList->front())));    

    for (const CartesianVector &splitPosition : {endpointProjection, vertexProjection})
    {
        ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
        std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

        // Only look at the most populated cluster
        const Cluster *const pCluster(clusterVector.front());
        this->SplitCluster(pCluster, splitPosition);
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionSimpleClusterSplittingAlgorithm::GetMCMuon(const MCParticleList *const pMCParticleList, const MCParticle *&pMCMuon) const
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if ((abs(pMCParticle->GetParticleId()) == 13) && (pMCParticle->GetParentList().empty()))
        {
            pMCMuon = pMCParticle;
            break; // Exit after finding the primary muon
        }
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionSimpleClusterSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, const CartesianVector &splitPosition) const
{
    CartesianVector splitAxis1(0.f, 0.f, 0.f), splitAxis2(0.f, 0.f, 0.f);
    if (!this->GetPathwayDirections(pCluster, splitPosition, splitAxis1, splitAxis2))
        return;

    CaloHitList caloHitList1, caloHitList2;
    this->SplitCaloHits(pCluster, splitPosition, splitAxis1, splitAxis2, caloHitList1, caloHitList2);

    if (caloHitList1.empty() || caloHitList2.empty())
        return;
    
    // Sanity check
    if ((LArClusterHelper::GetClosestDistance(splitPosition, caloHitList1) > m_maxDistToSplitPos) ||
        (LArClusterHelper::GetClosestDistance(splitPosition, caloHitList2) > m_maxDistToSplitPos))
    {
        std::cout << "Warning: Split cluster does not contain hits close to the muon endpoint!" << std::endl;
        return;
    }    

    // Split cluster
    const ClusterList clusterList(1, pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,PandoraContentApi::InitializeFragmentation(
        *this, clusterList, clusterListToDeleteName, clusterListToSaveName));

    for (const CaloHitList &caloHitList : {caloHitList1, caloHitList2})
    {
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList = caloHitList;

        const Cluster *pNewCluster(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
    }

    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SolutionSimpleClusterSplittingAlgorithm::GetPathwayDirections(const Cluster *const pCluster, const CartesianVector &splitPosition,
    CartesianVector &splitAxis1, CartesianVector &splitAxis2) const
{
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    // Bin hits wrt their angle relative to the x-axis
    std::map<int, int> angleToXAxis;    
    const float binInterval(2.f * M_PI / m_nAngularBins);

    for (const CaloHit *const pHit : clusterHits)
    {
        const CartesianVector hitPosition(pHit->GetPositionVector());
        const CartesianVector relativePosition(hitPosition - splitPosition);

        if (relativePosition.GetMagnitude() < std::numeric_limits<float>::epsilon()) { continue; }
        if (relativePosition.GetMagnitude() > m_radiusForDirEstimate) { continue; }
        
        float angle(std::atan2(relativePosition.GetZ(), relativePosition.GetX()));
        if (angle < 0.f) { angle += 2.f * M_PI; }
        const int binIndex(std::floor(angle / binInterval));
        
        if (angleToXAxis.find(binIndex) == angleToXAxis.end())
            angleToXAxis[binIndex] = 1;
        else
            ++angleToXAxis[binIndex];
    }

    // Find two highest populated bins
    std::vector<std::pair<int, int>> sortedBins(angleToXAxis.begin(), angleToXAxis.end());
    std::sort(sortedBins.begin(), sortedBins.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.second > b.second;
    });

    if (sortedBins.size() < 2)
        return false;

    splitAxis1 = CartesianVector(std::cos(sortedBins[0].first * binInterval), 0.f, std::sin(sortedBins[0].first * binInterval));
    splitAxis2 = CartesianVector(std::cos(sortedBins[1].first * binInterval), 0.f, std::sin(sortedBins[1].first * binInterval));
    
    /////// Visualisation /////
    //const CartesianVector splitAxis1End(splitPosition + (splitAxis1 * 10.f));
    //const CartesianVector splitAxis2End(splitPosition + (splitAxis2 * 10.f));
    //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitPosition, "Split Position", BLACK, 2));
    //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &splitPosition, &splitAxis1End, "Split Axis 1", RED, 2, 1));
    //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &splitPosition, &splitAxis2End, "Split Axis 2", BLUE, 2, 1));
    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    
    return true;
}

    
//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionSimpleClusterSplittingAlgorithm::SplitCaloHits(const Cluster *const pCluster, const CartesianVector &splitPosition, const CartesianVector &splitAxis1,
    const CartesianVector &splitAxis2, CaloHitList &caloHitList1, CaloHitList &caloHitList2) const
{
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    for (const CaloHit *const pHit : clusterHits)
    {        
        const CartesianVector relativePosition(pHit->GetPositionVector() - splitPosition);        
        const float l1(relativePosition.GetDotProduct(splitAxis1)), l2(relativePosition.GetDotProduct(splitAxis2));

        if (((l1 > 0.f) && (l2 > 0.f)) || ((l1 < 0.f) && (l2 < 0.f)))
        {
            const float t1(relativePosition.GetCrossProduct(splitAxis1).GetMagnitude()), t2(relativePosition.GetCrossProduct(splitAxis2).GetMagnitude());
            (t1 < t2) ? caloHitList1.push_back(pHit) : caloHitList2.push_back(pHit);
        }
        else if (l1 > 0.f)
        {
            caloHitList1.push_back(pHit);
        }
        else if (l2 > 0.f)
        {
            caloHitList2.push_back(pHit);
        }
    }

    //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList1, "CaloHitList1", RED));
    //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList2, "CaloHitList2", BLUE));
    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionSimpleClusterSplittingAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NAngularBins", m_nAngularBins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RadiusForDirEstimate", m_radiusForDirEstimate));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxDistToSplitPos", m_maxDistToSplitPos));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

