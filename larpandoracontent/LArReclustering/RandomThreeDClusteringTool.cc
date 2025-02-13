/**
 *  @file   larpandoracontent/LArReclustering/RandomThreeDClusteringTool.cc
 *
 *  @brief  Implementation file for the AI-enchanced data-driven AGI super pro mlg reclustering
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArReclustering/RandomThreeDClusteringTool.h"

#include <random>

using namespace pandora;

namespace lar_content
{

RandomThreeDClusteringTool::RandomThreeDClusteringTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool RandomThreeDClusteringTool::Run(const CaloHitList &inputCaloHitList, std::vector<CaloHitList> &outputCaloHitListsVector)
{
    std::mt19937 gen(0);
    std::uniform_int_distribution<> distr(1, inputCaloHitList.size());
    const int splitIdx {distr(gen)};

    // Lists for hits that have indices to the left and right of the random number
    CaloHitList leftCaloHitList, rightCaloHitList;

    int idx {0};
    for (const CaloHit *const pCaloHit3D : inputCaloHitList)
    {
        if (idx++ < splitIdx)
            leftCaloHitList.push_back(pCaloHit3D);
        else
            rightCaloHitList.push_back(pCaloHit3D);
    }
    std::cout << leftCaloHitList.size() << " -- " << rightCaloHitList.size() << "\n";
    outputCaloHitListsVector.push_back(leftCaloHitList);
    outputCaloHitListsVector.push_back(rightCaloHitList);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RandomThreeDClusteringTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
