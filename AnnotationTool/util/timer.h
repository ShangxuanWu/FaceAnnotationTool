#pragma once

#ifdef _OPENMP
#	include <omp.h>
#else
#	include <time.h>
#endif

#include <string>
#include <iostream>
#include <math.h>

class Timer {
public:
	Timer()
		: startTime(0)
	{
		startTime = elapsed();
	}

	double elapsed() const {
#ifdef _OPENMP
		return omp_get_wtime() - startTime;
#else
		return (clock() - startTime) / CLOCKS_PER_SEC;
#endif
	}

private:
	double startTime;
};

inline std::ostream& operator<<(std::ostream& os, const Timer& timer) {
	double total = timer.elapsed();
	size_t hours = size_t(floor(total / 3600));
	if (hours)
		os << hours << "h ";
	total -= hours * 3600;
	size_t mins = size_t(floor(total / 60));

	if (mins || hours)
		os << mins << "m ";
	total -= 60 * mins;

	return os << round(total * 100) / 100. << " s";
}

class AutoTimer : public Timer {
public:
	explicit AutoTimer(const char* str = "Total time elapsed", std::ostream& outStream = std::cerr)
		: str(str)
		, outStream(outStream)
	{
	}

	AutoTimer(const std::string& str, std::ostream& outStream = std::cerr)
		: str2(str)
		, str(str2.c_str())
		, outStream(outStream)
	{
	}

	~AutoTimer() {
		outStream << *this << ": " << str << std::endl;
	}

	std::ostream& stream() {
		return outStream;
	}
	
private:
	std::string str2;
	const char* str;
	std::ostream& outStream;
};
