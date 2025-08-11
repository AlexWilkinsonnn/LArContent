/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning shower growing algorithm
 *
 *  $Log: $
 */

#include "larpandoradlcontent/LArShowerGrowing/DLShowerGrowingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLShowerGrowingAlgorithm::DLShowerGrowingAlgorithm() :
    m_deltaRayLengthThresholdSquared{std::map<HitType, float>{}},
    m_deltaRayParentWeightThreshold{0.f},
    m_trainingMode{false},
    m_trainingTreeName{""}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLShowerGrowingAlgorithm::~DLShowerGrowingAlgorithm()
{
    if (m_trainingMode)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName, m_trainingFileName, "UPDATE"));

        for (const HitType &view : { TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W })
        {
            const float pitch{LArGeometryHelper::GetWirePitch(this->GetPandora(), view)};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_view_data", "view", static_cast<int>(view)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_view_data", "pitch", pitch));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName + "_view_data"));
        }
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName + "_view_data", m_trainingFileName, "UPDATE"));

    }
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
    // NOTE Other features:
    // - hit distances to nearest X and Z gaps (XXX only the cathode gap seems to be in the 1x2x6 geometry?)

    const VertexList *pVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    const Vertex *const pVertex{pVertexList->front()};
    if (pVertex->GetVertexType() != VERTEX_3D)
    {
        return STATUS_CODE_INVALID_PARAMETER;
    }
    std::map<HitType, CartesianVector> viewToVtxPos = {
        { TPC_VIEW_U, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U) },
        { TPC_VIEW_V, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V) },
        { TPC_VIEW_W, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W) }};

    // Gather clusters
    ClusterList clusterList;
    for (std::string listName : m_clusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pClusterList));
        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }

    int currClusterID{0}, currMCID{0};
    std::map<const MCParticle* const, int> mcToID = { { nullptr, -1 } }, mcToPDG = { { nullptr, 0 } };
    std::vector<int> clusterView, clusterID;
    std::vector<int> mcID = { -1 }, mcPDG = { 0 };
    std::vector<int> hitClusterID;
    std::vector<float> hitXRelPos, hitZRelPos;
    std::vector<float> hitRRelPos, hitSinThetaRelPos, hitCosThetaRelPos;  // Polar coordinates
    std::vector<float> hitXWidth;
    std::vector<float> hitDistToXGap;
    std::vector<float> hitEnergy;
    std::vector<int> hitMCID;

    std::map<const MCParticle *const, const MCParticle *const> mcFoldTo;

    std::set<double> xGaps;
    for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
    {
        const LineGap *const pLineGap = dynamic_cast<const LineGap *>(pDetectorGap);
        if (pLineGap->GetLineGapType() == TPC_DRIFT_GAP)
        {
            xGaps.insert(static_cast<double>(pLineGap->GetLineStartX()));
            xGaps.insert(static_cast<double>(pLineGap->GetLineEndX()));
        }
    }

    for (const Cluster *const pCluster : clusterList)
    {
        const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        clusterView.emplace_back(static_cast<int>(view));
        clusterID.emplace_back(currClusterID);

        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        CaloHitList clusterCaloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
        clusterCaloHits.insert(clusterCaloHits.end(), isolatedHits.begin(), isolatedHits.end());
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            const double x{static_cast<double>(pCaloHit->GetPositionVector().GetX())};
            const double xRel{x - static_cast<double>(viewToVtxPos.at(view).GetX())};
            const double z{static_cast<double>(pCaloHit->GetPositionVector().GetZ())};
            const double zRel{z - static_cast<double>(viewToVtxPos.at(view).GetZ())};

            const double rRel{std::sqrt(pow(xRel, 2.) + pow(zRel, 2.))};
            const double cosThetaRel{ rRel != 0. ? xRel / rRel : 0. };
            const double sinThetaRel{ rRel != 0. ? zRel / rRel : 0. };

            double distToXGap{std::numeric_limits<double>::max()};
            for (const double xGap : xGaps)
            {
                const double dist{x - xGap};
                if (std::abs(dist) < std::abs(distToXGap))
                {
                    distToXGap = dist;
                }
            }

            const MCParticle *const pMainMC{this->GetMainMC(pCaloHit, mcFoldTo)};
            if (mcToID.find(pMainMC) == mcToID.end())
            {
                mcToID.insert({ pMainMC, currMCID++ });
                mcToPDG.insert({ pMainMC, pMainMC->GetParticleId() });
                mcID.emplace_back(mcToID.at(pMainMC));
                mcPDG.emplace_back(mcToPDG.at(pMainMC));
            }

            hitClusterID.emplace_back(currClusterID);

            hitXRelPos.emplace_back(static_cast<float>(xRel));
            hitZRelPos.emplace_back(static_cast<float>(zRel));

            hitRRelPos.emplace_back(static_cast<float>(rRel));
            hitCosThetaRelPos.emplace_back(static_cast<float>(cosThetaRel));
            hitSinThetaRelPos.emplace_back(static_cast<float>(sinThetaRel));

            hitXWidth.emplace_back(pCaloHit->GetCellSize1());

            hitDistToXGap.emplace_back(static_cast<float>(distToXGap));

            hitEnergy.emplace_back(pCaloHit->GetMipEquivalentEnergy());

            hitMCID.emplace_back(mcToID.at(pMainMC));
        }
        currClusterID++;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_id", &clusterID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_view", &clusterView));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_id", &mcID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_pdg", &mcPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_cluster_id", &hitClusterID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_rel_pos", &hitXRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_z_rel_pos", &hitZRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_r_rel_pos", &hitRRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_ctheta_rel_pos", &hitCosThetaRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_stheta_rel_pos", &hitSinThetaRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_width", &hitXWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_gap_dist", &hitDistToXGap));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_energy", &hitEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_mc_id", &hitMCID));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* DLShowerGrowingAlgorithm::GetMainMC(
    const CaloHit *const pCaloHit, std::map<const MCParticle *const, const MCParticle *const> &mcFoldTo) const
{
    MCParticleWeightMap weightMap{pCaloHit->GetMCParticleWeightMap()};
    MCParticleWeightMap foldedWeightMap;
    for (const auto &[pMC, weight] : weightMap)
    {
        const MCParticle *pFoldedMC{nullptr};
        if (mcFoldTo.find(pMC) != mcFoldTo.end())
        {
            pFoldedMC = mcFoldTo.at(pMC);
        }
        else
        {
            pFoldedMC = this->FoldMCTo(pMC);
            mcFoldTo.insert({pMC, pFoldedMC});
        }
        foldedWeightMap[pFoldedMC] += weight;
    }
    weightMap = foldedWeightMap;

    const MCParticle *pMainMC{nullptr};
    float maxWeight{0.f};
    for (const auto &[pMC, weight] : weightMap)
    {
        if (weight > maxWeight)
        {
            pMainMC = pMC;
            maxWeight = weight;
        }
        if (weight == maxWeight) // tie-breaker (very unlikely)
        {
            if (LArMCParticleHelper::SortByMomentum(pMC, pMainMC))
            {
                pMainMC = pMC;
            }
        }
    }

    if (pMainMC)
    {
        pMainMC = this->FoldPotentialDeltaRayTo(pCaloHit, pMainMC);
    }

    return pMainMC;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* DLShowerGrowingAlgorithm::FoldMCTo(const MCParticle *const pMC) const
{
    if (!this->IsEM(pMC))
    {
        return pMC;
    }

    const MCParticle *pCurrentMC{pMC};
    const MCParticle *pLeadingMC{pMC};
    while (!pCurrentMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{*(pCurrentMC->GetParentList().begin())};
        const int parentPdg{std::abs(pParentMC->GetParticleId())};
        if (parentPdg == PHOTON || parentPdg == E_MINUS)
        {
            pCurrentMC = pParentMC;
            continue;
        }
        pLeadingMC = pCurrentMC;
        break;
    }

    return pLeadingMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* DLShowerGrowingAlgorithm::FoldPotentialDeltaRayTo(const CaloHit *const pCaloHit, const MCParticle *const pMC) const
{
    // Not an electron -> not a delta ray -> do nothing
    if (pMC->IsRootParticle() || pMC->GetParticleId() != E_MINUS)
    {
        return pMC;
    }

    // Did not come from a track-like particle -> not a delta ray -> do nothing
    const MCParticle *const pParentMC{*(pMC->GetParentList().begin())};
    const int parentPdg{std::abs(pParentMC->GetParticleId())};
    if (parentPdg == PHOTON || parentPdg == E_MINUS || PdgTable::GetParticleCharge(parentPdg) == 0)
    {
        return pMC;
    }

    // Delta ray that does not start a shower and is short -> fold into parent particle
    if (!this->CausesShower(pMC, 0) &&
        (pMC->GetVertex() - pMC->GetEndpoint()).GetMagnitudeSquared() < m_deltaRayLengthThresholdSquared.at(pCaloHit->GetHitType()))
    {
        return pParentMC;
    }

    // Now have a delta ray that we would like to cluster but only the hits that are not overlapping with the parent particle
    float parentWeight{std::numeric_limits<float>::lowest()};
    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
    for (const auto &[pContributingMC, weight] : weightMap)
    {
        if (pContributingMC == pParentMC)
        {
            parentWeight = weight;
            break;
        }
    }
    if (parentWeight > m_deltaRayParentWeightThreshold)
    {
        return pParentMC;
    }
    return pMC;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline bool DLShowerGrowingAlgorithm::IsEM(const pandora::MCParticle *const pMC) const
{
    const int pdg{std::abs(pMC->GetParticleId())};
    return (pdg == E_MINUS || pdg == PHOTON);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DLShowerGrowingAlgorithm::CausesShower(const MCParticle *const pMC, int nDescendentElectrons) const
{
    if (nDescendentElectrons > 1)
    {
        return true;
    }

    if (std::abs(pMC->GetParticleId()) == E_MINUS)
    {
        nDescendentElectrons++; // Including the parent particle, ie. the first in the recursion, as a descendent
    }
    for (const MCParticle *pChildMC : pMC->GetDaughterList())
    {
        if (this->CausesShower(pChildMC, nDescendentElectrons))
        {
            return true;
        }
    }

    return false;
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
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    m_deltaRayLengthThresholdSquared = {
        { TPC_VIEW_U, static_cast<float>(pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U), 2.)) },
        { TPC_VIEW_V, static_cast<float>(pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V), 2.)) },
        { TPC_VIEW_W, static_cast<float>(pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W), 2.)) }
    };

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
