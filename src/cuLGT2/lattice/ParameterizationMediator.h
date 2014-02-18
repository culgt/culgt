/**
 *
 *  Created on: Feb 11, 2014
 *      Author: vogt
 */

#ifndef PARAMETERIZATIONMEDIATOR_H_
#define PARAMETERIZATIONMEDIATOR_H_

namespace culgt
{

/*
 *
 */
template<typename ParamType1, typename ParamType2, typename LinkType1, typename LinkType2> class ParameterizationMediator
{
public:
	static void assign( LinkType1& l1, const LinkType2& l2 );
};

} /* namespace culgt */

#endif /* PARAMETERIZATIONMEDIATOR_H_ */
