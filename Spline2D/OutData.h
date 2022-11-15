#pragma once
#include <fstream>
#include <vector>
#include "programlog.h"
/// <summary>
/// Ñòðóêòóðà äëÿ âõîäíûõ äàííûõ
/// </summary>
struct OutData
{
	virtual std::ostream& writeData(std::ostream& out) = 0;
	virtual std::string getTypeDataName() = 0;
	std::ostream& operator >> (std::ostream& out) { writeData(out); }
};

//std::istream& operator >> (std::istream& in, InData& data);

template<class DataType, typename = std::enable_if_t<std::is_base_of<OutData, DataType>::value>>
std::ostream& operator << (std::ostream& in, DataType& data)
{
	return data.writeData(in);
}

template<class DataType, typename = std::enable_if_t<std::is_base_of<OutData, DataType>::value>>
std::ostream& writeData(std::ostream& out, std::vector<DataType>& data)
{
	DataType buf;

	int nData = data.size();
	out << nData << std::endl;

	if (out.fail())
	{
		programlog::writeErr("Cannot write number of \"" + buf.getTypeDataName() + "\" elements");
		return out;
	}


	for (int i = 0; i < nData; i++)
	{
		out << data[i] << std::endl;
		if (out.fail())
		{
			programlog::writeErr("Cannot write \"" + buf.getTypeDataName() + "\" element");
			return out;
		}
	}
	return out;
}

template<class DataType, typename = std::enable_if_t<std::is_base_of<OutData, DataType>::value>>
std::ostream& operator << (std::ostream& in, std::vector<DataType>& dataAr)
{
	return writeData(in, dataAr);
}

