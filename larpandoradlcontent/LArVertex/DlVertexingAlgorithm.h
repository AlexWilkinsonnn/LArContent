/**
 *  @file   larpandoradlcontent/LArVertexing/DlVertexingAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_VERTEXING_ALGORITHM_H
#define LAR_DL_VERTEXING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArObjects/VertexTuple.h"
#include "larpandoradlcontent/LArVertex/DlVertexingBaseAlgorithm.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningTrackShowerIdAlgorithm class
 */
class DlVertexingAlgorithm : public DlVertexingBaseAlgorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DlVertexingAlgorithm();

    ~DlVertexingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /*
     *  @brief  Create input for the network from a calo hit list
     *
     *  @param  caloHits The CaloHitList from which the input should be made
     *  @param  view The wire plane view
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     *  @param  networkInput The TorchInput object to populate
     *  @param  pixelVector The output vector of populated pixels
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeNetworkInputFromHits(const pandora::CaloHitList &caloHits, const pandora::HitType view, const float xMin,
        const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const;

    /*
     *  @brief  Create a list of wire plane-space coordinates from a canvas
     *
     *  @param  canvas The input canvas
     *  @param  canvasWidth The width of the canvas
     *  @param  canvasHeight The height of the canvas
     *  @param  columnOffset The column offset used when populating the canvas
     *  @param  rowOffset The row offset used when populating the canvas
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     *  @param  positionVector The output vector of wire plane positions
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeWirePlaneCoordinatesFromCanvas(float **canvas, const int canvasWidth, const int canvasHeight,
        const int columnOffset, const int rowOffset, const pandora::HitType view, const float xMin, const float xMax, const float zMin,
        const float zMax, pandora::CartesianPointVector &positionVector) const;

    /**
     *  @brief Create a vertex list from the candidate vertices.
     *
     *  @param  candidates The candidate positions with which to create the list.
     *
     *  @return The StatusCode resulting from the function
     */
    pandora::StatusCode MakeCandidateVertexList(const pandora::CartesianPointVector &positions);

#ifdef MONITORING
    /**
     *  @brief  Populate a root true with vertex information.
     */
    void PopulateRootTree(const std::vector<VertexTuple> &vertexTuples, const pandora::CartesianPointVector &vertexCandidatesU,
        const pandora::CartesianPointVector &vertexCandidatesV, const pandora::CartesianPointVector &vertexCandidatesW) const;
#endif

    int m_event;                ///< The current event number
    bool m_visualise;           ///< Whether or not to visualise the candidate vertices
    bool m_writeTree;           ///< Whether or not to write validation details to a ROOT tree
    std::string m_rootTreeName; ///< The ROOT tree name
    std::string m_rootFileName; ///< The ROOT file name
    std::mt19937 m_rng;         ///< The random number generator
};

} // namespace lar_dl_content

#endif // LAR_DL_VERTEXING_ALGORITHM_H
