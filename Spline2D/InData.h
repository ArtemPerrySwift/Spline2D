#pragma once
#include <fstream>
#include <vector>
#include "programlog.h"
/// <summary>
/// Ñòðóêòóðà äëÿ âõîäíûõ äàííûõ
/// </summary>
struct InData
{
	virtual std::istream& readData(std::istream& in) = 0;
	virtual std::string getTypeDataName() = 0;
	std::istream& operator << (std::istream& in) { readData(in); }
};

//std::istream& operator >> (std::istream& in, InData& data);

template<class DataType, typename = std::enable_if_t<std::is_base_of<InData, DataType>::value>>
std::istream& operator >> (std::istream& in, DataType& data)
{
	return data.readData(in);
}

template<class DataType, typename = std::enable_if_t<std::is_base_of<InData, DataType>::value>>
std::istream& readData(std::istream& in, std::vector<DataType>& data)
{
	DataType buf;

	int nData;
	in >> nData;

	if (in.fail())
	{
		programlog::writeErr("Cannot read number of \"" + buf.getTypeDataName() + "\" elements");
		return in;
	}

	data.resize(nData);
	if (nData < 0)
	{
		programlog::writeErr("The number of \"" + buf.getTypeDataName() + "\" elements cannot be < 0");
		return in;
	}


	for (int i = 0; i < nData; i++)
	{
		in >> buf;
		data[i] = buf;
		if (in.fail())
		{
			programlog::writeErr("Cannot read \"" + buf.getTypeDataName() + "\" element");
			return in;
		}
	}
	return in;
}

template<class DataType, typename = std::enable_if_t<std::is_base_of<InData, DataType>::value>>
std::istream& operator >> (std::istream& in, std::vector<DataType>& dataAr)
{
	return readData(in, dataAr);
}
