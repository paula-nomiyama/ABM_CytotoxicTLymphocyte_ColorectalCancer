#ifndef TCELLPROLIFERATIVETYPE_HPP_
#define TCELLPROLIFERATIVETYPE_HPP_

#include "AbstractCellProliferativeType.hpp"

class TCellProliferativeType : public AbstractCellProliferativeType
{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProliferativeType>(*this);
    }

public:
    TCellProliferativeType();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TCellProliferativeType)

#endif