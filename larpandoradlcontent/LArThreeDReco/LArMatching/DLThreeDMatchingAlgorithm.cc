/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLThreeDMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning shower growing algorithm
 *
 *  $Log: $
 */

#include <torch/torch.h>

#include "Geometry/LArReadoutVolume.h"
#include "Geometry/LArReadoutUnit.h"
#include "Geometry/LArReadoutChannel.h"

#include "larpandoradlcontent/LArThreeDReco/LArMatching/DLThreeDMatchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLThreeDMatchingAlgorithm::HitFeatures::HitFeatures() :
    m_xRel{0.f},
    m_zRel{0.f},
    m_rRel{0.f},
    m_cosThetaRel{0.f},
    m_sinThetaRel{0.f},
    m_distToXGap{0.f},
    m_xWidth{0.f},
    m_energy{0.f},
    m_uMinZRel{0.f},
    m_uMaxZRel{0.f},
    m_vMinZRel{0.f},
    m_vMaxZRel{0.f},
    m_wMinZRel{0.f},
    m_wMaxZRel{0.f},
    m_larTpcVolId{0},
    m_daughterVolId{0}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLThreeDMatchingAlgorithm::DLThreeDMatchingAlgorithm() :
    m_polarRScaleFactor{1.f},
    m_cartesianXScaleFactor{1.f},
    m_cartesianZScaleFactor{1.f},
    m_trainingMode{false},
    m_trainingTreeName{""}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLThreeDMatchingAlgorithm::~DLThreeDMatchingAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            this->WriteTrainingSample();
        }
        catch (StatusCodeException &e)
        {
            std::cout << "DLThreeDMatchingAlgorithm: Error while writing out training sample ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDMatchingAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDMatchingAlgorithm::PrepareTrainingSample() const
{
    m_rollUpper.Reset();

    const std::map<HitType, CartesianVector> viewToVtxPos{this->Get2DVertices()};

    const ClusterList clusterList{this->GetAllClusters()};

    int currClusterID{0}, currMCID{0};
    std::map<const MCParticle *const, int> mcToID = {{nullptr, -1}}; // Initialise for null MCParticle
    std::vector<int> clusterView, clusterID;
    std::vector<int> mcID = {-1}, mcPDG = {0}, mcIsFromBeam = {0}, mcProcess = {0}; // Initialise for null MCParticle
    std::vector<int> hitClusterID;
    std::vector<float> hitXRelPos, hitZRelPos;
    std::vector<float> hitRRelPos, hitSinThetaRelPos, hitCosThetaRelPos; // Polar coordinates
    std::vector<float> hitXWidth;
    std::vector<float> hitDistToXGap;
    std::vector<float> hitEnergy;
    std::vector<float> hitUMinZRel, hitUMaxZRel, hitVMinZRel, hitVMaxZRel, hitWMinZRel, hitWMaxZRel;
    std::vector<int> hitMCID;
    std::vector<float> hitMCWeight;
    std::vector<int> hitLArTPCVolID, hitDaughterVolID;

    for (const Cluster *const pCluster : clusterList)
    {
        clusterView.emplace_back(static_cast<int>(LArClusterHelper::GetClusterHitType(pCluster)));
        clusterID.emplace_back(currClusterID);

        CaloHitList clusterCaloHits;
        LArClusterHelper::GetAllHits(pCluster, clusterCaloHits);
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            HitFeatures hitFeatures;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CalculateHitFeatures(pCaloHit, viewToVtxPos, hitFeatures));

            const MCParticle *const pMainMC{m_rollUpper.RollUpCaloHit(pCaloHit)};
            if (mcToID.find(pMainMC) == mcToID.end())
            {
                mcToID.insert({pMainMC, currMCID++});
                mcID.emplace_back(mcToID.at(pMainMC));
                mcPDG.emplace_back(pMainMC->GetParticleId());

                const MCParticle *const pMainMCParent{LArMCParticleHelper::GetParentMCParticle(pMainMC)};
                mcIsFromBeam.emplace_back(
                    (LArMCParticleHelper::IsBeamParticle(pMainMCParent) || LArMCParticleHelper::IsNeutrino(pMainMCParent)) ? 1 : 0);

                const LArMCParticle *pLArMainMC{dynamic_cast<const LArMCParticle *>(pMainMC)};
                PANDORA_RETURN_IF(STATUS_CODE_NOT_FOUND, !pLArMainMC);
                mcProcess.emplace_back(pLArMainMC->GetProcess());
            }

            hitClusterID.emplace_back(currClusterID);

            hitXRelPos.emplace_back(hitFeatures.m_xRel);
            hitZRelPos.emplace_back(hitFeatures.m_zRel);
            hitRRelPos.emplace_back(hitFeatures.m_rRel);
            hitCosThetaRelPos.emplace_back(hitFeatures.m_cosThetaRel);
            hitSinThetaRelPos.emplace_back(hitFeatures.m_sinThetaRel);
            hitXWidth.emplace_back(hitFeatures.m_xWidth);
            hitDistToXGap.emplace_back(hitFeatures.m_distToXGap);
            hitEnergy.emplace_back(hitFeatures.m_energy);
            hitUMinZRel.emplace_back(hitFeatures.m_uMinZRel);
            hitUMaxZRel.emplace_back(hitFeatures.m_uMaxZRel);
            hitVMinZRel.emplace_back(hitFeatures.m_vMinZRel);
            hitVMaxZRel.emplace_back(hitFeatures.m_vMaxZRel);
            hitWMinZRel.emplace_back(hitFeatures.m_wMinZRel);
            hitWMaxZRel.emplace_back(hitFeatures.m_wMaxZRel);
            hitLArTPCVolID.emplace_back(static_cast<int>(hitFeatures.m_larTpcVolId));
            hitDaughterVolID.emplace_back(static_cast<int>(hitFeatures.m_daughterVolId));

            hitMCID.emplace_back(mcToID.at(pMainMC));

            const MCParticleWeightMap weightMap{m_rollUpper.GetRolledUpMCParticleWeightMap(pCaloHit)};
            PANDORA_RETURN_IF(STATUS_CODE_NOT_ALLOWED, weightMap.empty() && pMainMC);
            hitMCWeight.emplace_back(pMainMC ? weightMap.at(pMainMC) : 1.f);
        }
        currClusterID++;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_id", &clusterID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "cluster_view", &clusterView));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_id", &mcID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_pdg", &mcPDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_is_from_beam", &mcIsFromBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "mc_process", &mcProcess));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_cluster_id", &hitClusterID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_rel_pos", &hitXRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_z_rel_pos", &hitZRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_r_rel_pos", &hitRRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_ctheta_rel_pos", &hitCosThetaRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_stheta_rel_pos", &hitSinThetaRelPos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_width", &hitXWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_x_gap_dist", &hitDistToXGap));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_energy", &hitEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_u_min_z_rel_pos", &hitUMinZRel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_u_max_z_rel_pos", &hitUMaxZRel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_v_min_z_rel_pos", &hitVMinZRel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_v_max_z_rel_pos", &hitVMaxZRel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_w_min_z_rel_pos", &hitWMinZRel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_w_max_z_rel_pos", &hitWMaxZRel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_mc_id", &hitMCID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_mc_weight", &hitMCWeight));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_tpc_vol_id", &hitLArTPCVolID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName, "hit_daughter_vol_id", &hitDaughterVolID));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DLThreeDMatchingAlgorithm::WriteTrainingSample() const
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName, m_trainingFileName, "UPDATE"));

    for ([[maybe_unused]] const HitType &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_view_data", "view", static_cast<int>(view)));
        PANDORA_MONITORING_API(SetTreeVariable(
            this->GetPandora(), m_trainingTreeName + "_view_data", "pitch", LArGeometryHelper::GetWirePitch(this->GetPandora(), view)));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName + "_view_data"));
    }
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName + "_view_data", m_trainingFileName, "UPDATE"));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_event_data", "scalefactor_polar_r", m_polarRScaleFactor))
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_event_data", "scalefactor_cartesian_x", m_cartesianXScaleFactor))
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_trainingTreeName + "_event_data", "scalefactor_cartesian_z", m_cartesianZScaleFactor))
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_trainingTreeName + "_event_data"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_trainingTreeName + "_event_data", m_trainingFileName, "UPDATE"));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::map<HitType, CartesianVector> DLThreeDMatchingAlgorithm::Get2DVertices() const
{
    const VertexList *pVertexList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    const Vertex *const pVertex{pVertexList->front()};
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, pVertex->GetVertexType() != VERTEX_3D);
    return {{TPC_VIEW_U, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U)},
        {TPC_VIEW_V, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V)},
        {TPC_VIEW_W, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W)}};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

ClusterList DLThreeDMatchingAlgorithm::GetAllClusters() const
{
    ClusterList clusterList;
    for (const std::string &listName : m_clusterListNames)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetClusters(listName, clusterList));
    }
    return clusterList;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDMatchingAlgorithm::GetClusters(const std::string clusterListName, ClusterList &clusterList) const
{
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));
    clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDMatchingAlgorithm::CalculateHitFeatures(
    const CaloHit *const pCaloHit, const std::map<HitType, CartesianVector> &viewToVtxPos, HitFeatures &hitFeatures) const
{
    const HitType view{pCaloHit->GetHitType()};
    const double x{static_cast<double>(pCaloHit->GetPositionVector().GetX())};
    const double xRel{x - static_cast<double>(viewToVtxPos.at(view).GetX())};
    const double z{static_cast<double>(pCaloHit->GetPositionVector().GetZ())};
    const double zRel{z - static_cast<double>(viewToVtxPos.at(view).GetZ())};

    const double rRel{std::sqrt(pow(xRel, 2.) + pow(zRel, 2.))};
    const double cosThetaRel{rRel != 0. ? xRel / rRel : 0.};
    const double sinThetaRel{rRel != 0. ? zRel / rRel : 0.};

    double distToXGap{std::numeric_limits<double>::max()};
    for (const double xGap : m_detectorXGaps)
    {
        const double dist{x - xGap};
        if (std::abs(dist) < std::abs(distToXGap))
            distToXGap = dist;
    }

    hitFeatures.m_xRel = static_cast<float>(xRel);
    hitFeatures.m_zRel = static_cast<float>(zRel);
    hitFeatures.m_rRel = static_cast<float>(rRel);
    hitFeatures.m_cosThetaRel = static_cast<float>(cosThetaRel);
    hitFeatures.m_sinThetaRel = static_cast<float>(sinThetaRel);
    hitFeatures.m_distToXGap = static_cast<float>(distToXGap);
    hitFeatures.m_xWidth = pCaloHit->GetCellSize1();
    hitFeatures.m_energy = pCaloHit->GetMipEquivalentEnergy();

    const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));
    PANDORA_RETURN_IF(STATUS_CODE_FAILURE, !pLArCaloHit);
    hitFeatures.m_larTpcVolId = pLArCaloHit->GetLArTPCVolumeId();
    hitFeatures.m_daughterVolId = pLArCaloHit->GetDaughterVolumeId();

    const LArReadoutVolume &readoutVol{LArGeometryHelper::GetReadoutVolume(this->GetPandora(), *pLArCaloHit)};
    const LArReadoutChannel &readoutCh{LArGeometryHelper::GetReadoutChannel(this->GetPandora(), *pLArCaloHit)};
    
    for (const LArReadoutUnit &readoutUnit : readoutVol.GetReadoutUnits())
    {
        // The channel intervals only cover the other views, in the hit's own view the interval is just the hit's channel
        const HitType unitView{readoutUnit.GetView()};
        const LArReadoutChannel::ChannelInterval chInterval{unitView == view
                ? LArReadoutChannel::ChannelInterval{pLArCaloHit->GetChannelId(), pLArCaloHit->GetChannelId()}
                : readoutCh.GetChannelInterval(unitView)};
        const unsigned int chLo{chInterval.first}, chHi{chInterval.second};
        const float zLo{readoutUnit.GetReferenceCoordinate() + readoutUnit.GetPitch() * (static_cast<float>(chLo) - 0.5f)};
        const float zHi{readoutUnit.GetReferenceCoordinate() + readoutUnit.GetPitch() * (static_cast<float>(chHi) + 0.5f)};

        // The pitch is signed, so zLo is not necessarily the lower bound
        const float vtxZ{viewToVtxPos.at(unitView).GetZ()};
        const float minZRel{std::min(zLo, zHi) - vtxZ};
        const float maxZRel{std::max(zLo, zHi) - vtxZ};

        PANDORA_RETURN_IF(STATUS_CODE_OUT_OF_RANGE, unitView == view && (zRel < minZRel || zRel > maxZRel)); // Sanity check

        switch (unitView)
        {
            case TPC_VIEW_U:
                hitFeatures.m_uMinZRel = minZRel;
                hitFeatures.m_uMaxZRel = maxZRel;
                break;
            case TPC_VIEW_V:
                hitFeatures.m_vMinZRel = minZRel;
                hitFeatures.m_vMaxZRel = maxZRel;
                break;
            case TPC_VIEW_W:
                hitFeatures.m_wMinZRel = minZRel;
                hitFeatures.m_wMaxZRel = maxZRel;
                break;
            default:
                return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDMatchingAlgorithm::Infer()
{
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeDMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingTreeName", m_trainingTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
        const std::map<HitType, float> lengthThresholds{{TPC_VIEW_U, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)},
            {TPC_VIEW_V, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)},
            {TPC_VIEW_W, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)}};
        m_rollUpper = RollUpper(std::make_unique<RollUpEMWithComptonFilterAndAmbiguousDeltaRayHitsPolicy>(0.f, lengthThresholds));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    if (m_trainingMode)
    {
        const LArGeometryHelper::DetectorBoundaries detBounds{LArGeometryHelper::GetDetectorBoundaries(this->GetPandora())};
        double xLow{static_cast<double>(detBounds.m_xBoundaries.first)}, xHigh{static_cast<double>(detBounds.m_xBoundaries.second)};
        // double yLow{static_cast<double>(detBounds.m_yBoundaries.first)}, yHigh{static_cast<double>(detBounds.m_yBoundaries.second)};
        double zLow{static_cast<double>(detBounds.m_zBoundaries.first)}, zHigh{static_cast<double>(detBounds.m_zBoundaries.second)};
        m_polarRScaleFactor = static_cast<float>(1. / std::sqrt(std::pow(xHigh - xLow, 2.) + std::pow(zHigh - zLow, 2.)));
        m_cartesianXScaleFactor = static_cast<float>(1. / (xHigh - xLow));
        m_cartesianZScaleFactor = static_cast<float>(1. / (zHigh - zLow));
        std::cout << "Scalefactors for this geometry\n"
                  << "Polar R:     " << m_polarRScaleFactor << "\n"
                  << "Cartesian X: " << m_cartesianXScaleFactor << "\n"
                  << "Cartesian Z: " << m_cartesianZScaleFactor << "\n";
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PolarRScaleFactor", m_polarRScaleFactor));
        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CartesianXScaleFactor", m_cartesianXScaleFactor));
        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CartesianZScaleFactor", m_cartesianZScaleFactor));
    }

    for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
    {
        const LineGap *const pLineGap = dynamic_cast<const LineGap *>(pDetectorGap);
        if (pLineGap && pLineGap->GetLineGapType() == TPC_DRIFT_GAP)
        {
            m_detectorXGaps.insert(static_cast<double>(pLineGap->GetLineStartX()));
            m_detectorXGaps.insert(static_cast<double>(pLineGap->GetLineEndX()));
        }
    }
    // PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, !m_trainingMode && m_includeDistToXGapFeature && m_detectorXGaps.empty());

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
