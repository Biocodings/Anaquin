#ifndef GI_RESOURCES_HPP
#define GI_RESOURCES_HPP

#include <string>

namespace Spike
{
    struct ChromosomeInfo
    {
        float v;

        // Eg: chrT
        std::string id;
    };
    
    struct MixtureInfo
    {
        float va;
        float vb;
    };

    struct Resources
    {
        static MixtureInfo mixture();
        static ChromosomeInfo chromo();
    };
}

#endif