/**
 *  @file   larpandoracontent/LArReclustering/RandomThreeDClusteringTool.h
 *
 *  @brief  Implementation file for the AI-enchanced data-driven AGI super pro mlg reclustering
 *
 *  $Log: $
 */

#ifndef LAR_RANDOM_THREE_D_CLUSTERING_TOOL_H
#define LAR_RANDOM_THREE_D_CLUSTERING_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

namespace lar_content
{

class RandomThreeDClusteringTool : public ClusteringTool
{
public:
    RandomThreeDClusteringTool();

private:
    bool Run(const pandora::CaloHitList &inputCaloHitList, std::vector<pandora::CaloHitList> &outputCaloHitListsVector);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/);
};

} // namespace lar_content

#endif // #endif LAR_RANDOM_THREE_D_CLUSTERING_TOOL_H
