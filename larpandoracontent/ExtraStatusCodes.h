#include <exception>
#include <string>

#define PANDORA_RETURN(StatusCode)                                                                      \
{                                                                                                       \
    std::cout << "    in function: " << __FUNCTION__ << std::endl;                                      \
    std::cout << "    in file:     " << __FILE__ << " line#: " << __LINE__ << std::endl;                \
    return StatusCode;                                                                                  \
}

#define PANDORA_RETURN_IF(StatusCode, Condition)                                                        \
{                                                                                                       \
    if (Condition)                                                                                      \
    {                                                                                                   \
        std::cout << #Condition << " return " << StatusCodeToString(StatusCode) << std::endl;           \
        std::cout << "    in function: " << __FUNCTION__ << std::endl;                                  \
        std::cout << "    in file:     " << __FILE__ << " line#: " << __LINE__ << std::endl;            \
        return StatusCode;                                                                              \
    }                                                                                                   \
}

#define PANDORA_THROW(StatusCode)                                                                       \
{                                                                                                       \
    std::cout << "    in function: " << __FUNCTION__ << std::endl;                                      \
    std::cout << "    in file:     " << __FILE__ << " line#: " << __LINE__ << std::endl;                \
    throw pandora::StatusCodeException(StatusCode);                                                     \
}

#define PANDORA_THROW_IF(StatusCode, Condition)                                                         \
{                                                                                                       \
    if (Condition)                                                                                      \
    {                                                                                                   \
        std::cout << #Condition << " throw " << StatusCodeToString(StatusCode) << std::endl;            \
        std::cout << "    in function: " << __FUNCTION__ << std::endl;                                  \
        std::cout << "    in file:     " << __FILE__ << " line#: " << __LINE__ << std::endl;            \
        throw pandora::StatusCodeException(StatusCode);                                                 \
    }                                                                                                   \
}
