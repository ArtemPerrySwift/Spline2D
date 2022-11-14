#pragma once
#include"indata.h"
#define USER_COORDS true

struct Coord2D : InData
{
	double x;
	double y;

	Coord2D() { x = y = 0; }

	Coord2D(double x, double y)
	{
		this->x = x;
		this->y = y;
	}

	void transIntoMas(double mas[2])
	{
		mas[0] = x;
		mas[1] = y;
	}

	std::istream& readData(std::istream& in) override
	{
		in >> x >> y;
		if (in.fail())
		{
			programlog::writeErr("Unable to read " + getTypeDataName());
		}
		return in;
	}

	std::string getTypeDataName() override
	{
		return "CoordXY";
	}

	Coord2D operator +(const Coord2D& coord2)
	{
		Coord2D newCoord;
		newCoord.x = this->x + coord2.x;
		newCoord.y = this->y + coord2.y;
		return newCoord;
	}

	Coord2D operator -(const Coord2D& coord2)
	{
		Coord2D newCoord;
		newCoord.x = this->x - coord2.x;
		newCoord.y = this->y - coord2.y;
		return newCoord;
	}

	Coord2D operator /(const Coord2D& coord2)
	{
		Coord2D newCoord;
		newCoord.x = this->x / coord2.x;
		newCoord.y = this->y / coord2.y;
		return newCoord;
	}

	Coord2D operator *(const Coord2D& coord2)
	{
		Coord2D newCoord;
		newCoord.x = this->x * coord2.x;
		newCoord.y = this->y * coord2.y;
		return newCoord;
	}

	Coord2D operator /(const double& a)
	{
		Coord2D newCoord;
		newCoord.x = this->x / a;
		newCoord.y = this->y / a;
		return newCoord;
	}

	Coord2D operator *(const double& a)
	{
		Coord2D newCoord;
		newCoord.x = this->x * a;
		newCoord.y = this->y * a;
		return newCoord;
	}

	Coord2D& operator +=(const Coord2D& coord2)
	{
		Coord2D newCoord;
		this->x += coord2.x;
		this->y += coord2.y;
		return *this;
	}

	Coord2D& operator *=(const Coord2D& coord2)
	{
		Coord2D newCoord;
		this->x *= coord2.x;
		this->y *= coord2.y;
		return *this;
	}

	Coord2D& operator *=(const double& a)
	{
		Coord2D newCoord;
		this->x *= a;
		this->y *= a;
		return *this;
	}

	Coord2D& operator =(const double a[2])
	{
		Coord2D newCoord;
		this->x = a[0];
		this->y = a[1];
		return *this;
	}
}; #pragma once
