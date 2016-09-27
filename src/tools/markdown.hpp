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
        
            typedef std::string Code;
        
            struct Item
            {
                virtual Code generate() const = 0;
            };
        
            class RItem : public Item
            {
                public:
                    RItem(const Code &code) : _code(code) {}
            
                    virtual Code generate() const
                    {
                        throw "";
                    }
            
                private:
                    const Code _code;
            };
        
            class TextItem : public Item
            {
                public:
                    TextItem(const Title &title, const Code &txt) : _txt(txt), _title(title) {}
            
                    virtual Scripts generate() const;

                private:
                    const Scripts _txt;
                    const Title   _title;
            };
        
            class Section
            {
                public:
                    Section(const Title &title) : _title(title) {}

                    inline void addText(const Title &title, const Scripts &txt)
                    {
                        _items.push_back(std::shared_ptr<Item>(new TextItem(title, txt)));
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
