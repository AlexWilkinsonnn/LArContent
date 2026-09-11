/**
 *  @file   larpandoracontent/LArSolutions/SolutionExampleAlgorithm.h
 *
 *  @brief  Header file for the SolutionExampleAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SOLUTION_EXAMPLE_H
#define LAR_SOLUTION_EXAMPLE_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SolutionExampleAlgorithm class
 */
class SolutionExampleAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SolutionExampleAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_favouriteNumber;
    std::string m_mcParticleListName;
};

} // namespace lar_content

#endif // #ifndef LAR_SOLUTION_EXAMPLE_H

