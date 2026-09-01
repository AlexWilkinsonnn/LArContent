#!/bin/bash

set -e

export MY_TEST_AREA

if [ -z "${MY_TEST_AREA}" ]; then
    echo "Error: MY_TEST_AREA is not set, run setup script."
    exit 1
fi

ALGORITHM_BASE="${1}"
ALGORITHM_NAME="${ALGORITHM_BASE}Algorithm"
HEADER_THING=$(echo "${ALGORITHM_BASE}" | sed -E 's/([a-z0-9])([A-Z])/\1_\2/g' | tr '[:lower:]' '[:upper:]')
ALGORITHM_PATH="${2#./}"
ALGORITHM_PATH="${ALGORITHM_PATH%/}"
LOCAL_ALGORITHM_PATH="larpandoracontent/${ALGORITHM_PATH}"
GLOBAL_ALGORITHM_PATH="${MY_TEST_AREA}/LArContent/larpandoracontent/${ALGORITHM_PATH}"
LOCAL_ALGORITHM_PATH="${LOCAL_ALGORITHM_PATH%/}"
GLOBAL_ALGORITHM_PATH="${GLOBAL_ALGORITHM_PATH%/}"

LOCAL_SOURCE_FILE="${LOCAL_ALGORITHM_PATH}/${ALGORITHM_NAME}.cc"
GLOBAL_SOURCE_FILE="${GLOBAL_ALGORITHM_PATH}/${ALGORITHM_NAME}.cc"
LOCAL_HEADER_FILE="${LOCAL_ALGORITHM_PATH}/${ALGORITHM_NAME}.h"
GLOBAL_HEADER_FILE="${GLOBAL_ALGORITHM_PATH}/${ALGORITHM_NAME}.h"
CMAKE_FILE="${MY_TEST_AREA}/LArContent/cmake/LArContent_sources.cmake"

echo "Algorithm name: ${ALGORITHM_NAME}"
echo "Local algorithm path: ${LOCAL_ALGORITHM_PATH}"
echo "Global algorithm path: ${GLOBAL_ALGORITHM_PATH}"
echo "LOCAL_SOURCE_FILE: ${LOCAL_SOURCE_FILE}" 
echo "GLOBAL_SOURCE_FILE: ${GLOBAL_SOURCE_FILE}" 

# Check that the required arguments were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <algorithm_name> <path>"
    exit 1
fi

# Create the directory if it doesn't already exist
mkdir -p "${GLOBAL_ALGORITHM_PATH}"

# Don't overwrite an existing source file
if [ -e "${GLOBAL_SOURCE_FILE}" ]; then
    echo "Error: ${GLOBAL_SOURCE_FILE} already exists."
    exit 1
fi

# Don't overwrite an existing header file
if [ -e "${GLOBAL_HEADER_FILE}" ]; then
    echo "Error: ${GLOBAL_HEADER_FILE} already exists."
    exit 1
fi

# ---------------------------------------------------------------------------
# .h
# ---------------------------------------------------------------------------

cat > ${GLOBAL_HEADER_FILE} <<EOF
/**
 *  @file   ${LOCAL_HEADER_FILE}
 *
 *  @brief  Header file for the ${ALGORITHM_NAME} class.
 *
 *  \$Log: \$
 */
#ifndef LAR_${HEADER_THING}_H
#define LAR_${HEADER_THING}_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ${ALGORITHM_NAME} class
 */
class ${ALGORITHM_NAME} : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ${ALGORITHM_NAME}();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_${HEADER_THING}_H

EOF

# ---------------------------------------------------------------------------
# .cc
# ---------------------------------------------------------------------------

cat > ${GLOBAL_SOURCE_FILE} <<EOF
/**
 *  @file   ${LOCAL_SOURCE_FILE}
 *
 *  @brief  Implementation of the ${ALGORITHM_NAME} class.
 *
 *  \$Log: \$
 */

#include "Pandora/AlgorithmHeaders.h"

#include "${LOCAL_HEADER_FILE}"

using namespace pandora;

namespace lar_content
{

${ALGORITHM_NAME}::${ALGORITHM_NAME}()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ${ALGORITHM_NAME}::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ${ALGORITHM_NAME}::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

EOF

# ---------------------------------------------------------------------------
# Add the source file to LAR_CONTENT_SRCS and keep it alphabetically sorted
# ---------------------------------------------------------------------------

python3 - <<EOF
from pathlib import Path

cmake_file = Path("${CMAKE_FILE}")
new_source = "    ${LOCAL_SOURCE_FILE}"

text = cmake_file.read_text()

start_marker = "set(LAR_CONTENT_SRCS\n"
start = text.index(start_marker) + len(start_marker)
end = text.index("\n)", start)

sources = text[start:end].splitlines()

# Only add the source if it isn't already present
if new_source not in sources:
    sources.append(new_source)
    sources.sort()

    text = text[:start] + "\n".join(sources) + text[end:]
    cmake_file.write_text(text)

    print(f"Added {new_source} to {cmake_file}")
else:
    print(f"{new_source} is already present in {cmake_file}")

EOF

# ---------------------------------------------------------------------------
# Now add algorithm to LArContent.cc
# ---------------------------------------------------------------------------

LAR_CONTENT_FILE="${MY_TEST_AREA}/LArContent/larpandoracontent/LArContent.cc"

python3 - << EOF

import sys

with open('${LAR_CONTENT_FILE}') as f:
    lines = f.readlines()

include_line = f'#include "${LOCAL_HEADER_FILE}"\n'

if include_line not in lines:
    # Find all larpandoracontent includes
    include_indices = [
        i for i, line in enumerate(lines)
        if line.startswith('#include "larpandoracontent/')
    ]

    # Find the first include alphabetically after ours
    inserted = False

    for i in include_indices:
        existing = lines[i].strip()

        if existing > include_line.strip():
            lines.insert(i, include_line)
            inserted = True
            break

    # If ours is alphabetically last, put it after the final include
    if not inserted:
        lines.insert(include_indices[-1] + 1, include_line)


# Add algorithm to LAR_ALGORITHM_LIST in alphabetical order
algorithm_registration = f"LAr${1}"

NAME_COLUMN = 45
CONTINUATION_COLUMN = 128

prefix = f'    d("{algorithm_registration}",'
new_line = prefix.ljust(NAME_COLUMN)
new_line += "${ALGORITHM_NAME})"
new_line = new_line.ljust(CONTINUATION_COLUMN)
new_line += chr(92) + chr(10)

if not any(f'd("LAr${1}",' in line for line in lines):

    start = next(
        i for i, line in enumerate(lines)
        if line.startswith("#define LAR_ALGORITHM_LIST")
    )

    end = next(
        i for i in range(start + 1, len(lines))
        if lines[i].startswith("#define LAR_ALGORITHM_TOOL_LIST")
    )

    inserted = False

    for i in range(start + 1, end):
        line = lines[i]

        if 'd("' not in line:
            continue

        existing_name = line.split('d("', 1)[1].split('"', 1)[0]

        if "LAr${1}" < existing_name:
            lines.insert(i, new_line)
            inserted = True
            break

    if not inserted:
        lines.insert(end, new_line)

with open('${LAR_CONTENT_FILE}', "w") as f:
    f.writelines(lines)
EOF
