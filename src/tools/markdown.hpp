#ifndef MARKDOWN_HPP
#define MARKDOWN_HPP

#include <string>
#include <vector>
#include "data/data.hpp"

namespace Anaquin
{
    typedef std::string Title;
    
    class MarkDown
    {
        private:
        
            class Item
            {
                public:
                    Item(const Title &title, const Scripts &txt) : _title(title), _txt(txt) {}
                
                    virtual Scripts generate() const = 0;
                
                protected:

                    // Eg: Gene Expression
                    const Title _title;
                
                    // Eg: R-script
                    const Scripts _txt;
            };
        
            struct RItem : public Item
            {
                RItem(const Title &title, const Scripts &txt) : Item(title, txt) {}
                
                virtual Scripts generate() const;
            };
        
            struct TextItem : public Item
            {
                TextItem(const Title &title, const Scripts &txt) : Item(title, txt) {}

                virtual Scripts generate() const;
            };
        
            class Section
            {
                public:
                    Section(const Title &title) : _title(title) {}

                    inline void addText(const Title &title, const Scripts &txt)
                    {
                        _items.push_back(std::shared_ptr<Item>(new TextItem(title, txt)));
                    }
            
                    inline void addRCode(const Title &title, const Scripts &txt)
                    {
                        _items.push_back(std::shared_ptr<Item>(new RItem(title, txt)));
                    }

                    inline Scripts generate() const;

                private:
                    const Title _title;            
                    std::vector<std::shared_ptr<Item>> _items;
            };
        
        public:

            inline void start(const Title &title)
            {
                _current = std::shared_ptr<Section>(new Section(title));
            }
        
            inline void end()
            {
                _sections.push_back(*_current);
                _current = nullptr;
            }

            inline void addRCode(const Title &title, const Scripts &txt)
            {
                _current->addRCode(title, txt);
            }

            inline void addText(const Title &title, const Scripts &txt)
            {
                _current->addText(title, txt);
            }

            Scripts generate(const Title &title) const;

        private:
        
            // Current section be modified
            std::shared_ptr<Section> _current;
        
            // Completed sections
            std::vector<Section> _sections;
    };
}

#endif
