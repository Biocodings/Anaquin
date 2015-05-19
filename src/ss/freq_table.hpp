#ifndef QQ_FREQUENCY_TABLE_HPP
#define QQ_FREQUENCY_TABLE_HPP

#include <map>
#include <algorithm>
#include "QData/Variables/Code.hpp"

namespace QQ
{
	template <typename T> struct FrequencyTable
	{
		FrequencyTable(const CodeFrame &frame, const Values &values);

		SampleSize total;

		T min, max;

		std::map<T, Count> counts;
		std::map<T, Percentage> percents;

		Real overallPercents;
	};

	template<typename T> FrequencyTable<T>::FrequencyTable(const CodeFrame &frame, const Values &values)
	{
		typedef CodeFrame::value_type Type;

		if (values.empty())
		{
			throw std::runtime_error("Failed to initialize a frequency table from an empty input");
		}

		total = 0;
		//min = std::numeric_limits<T>::max();
		//max = std::numeric_limits<T>::min();

		std::for_each(frame.begin(), frame.end(), [&](Code i)
		{
			counts[i.value()] = 0;
		});

		std::for_each(values.begin(), values.end(), [&](CodeValue i)
		{
			total++;
			counts[i]++;
			//min = std::min(i, min);
			//max = std::max(i, max);
		});

		overallPercents = 0;

		std::for_each(counts.begin(), counts.end(), [&](std::pair<T, Count> p)
		{
			percents[p.first] = (1.0 * p.second) / total;
			overallPercents += percents[p.first];
		});

		overallPercents = overallPercents / counts.size();
	}
}

#endif