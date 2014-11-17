/*
 * SequencerRunner.h
 *
 *  Created on: Nov 7, 2014
 *      Author: vogt
 */

#ifndef SEQUENCERUNNER_H_
#define SEQUENCERUNNER_H_


#include <boost/mpl/if.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/unpack_args.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/lambda.hpp>

namespace mpl = boost::mpl;
using namespace std;
using namespace mpl::placeholders;

struct NIL
{
public:
	static const int value = 0;
};


//struct printseq
//{
//	template<typename T> void operator()(T)
//	{
//		cout << typeid(T).name() << "/";
//	}
//};



template<typename Chooser, typename Seq, typename T, typename... Ts> class SequenceRunner
{
public:
	template<typename S> void operator()(S)
	{
		typedef typename mpl::if_< mpl::bool_<(sizeof...(Ts) > 1)>, SequenceRunner<Chooser, typename mpl::push_back<Seq,S>::type,Ts...>, SequenceRunner<Chooser, typename mpl::push_back<Seq,S>::type,Ts...,NIL> >::type VSub;
		mpl::for_each<T>( VSub() );
	}
};

template<typename Chooser, typename Seq, typename T> class SequenceRunner<Chooser, Seq, T, NIL>
{
public:
	template<typename S> void operator()(S)
	{
		typedef SequenceRunner<Chooser, typename mpl::push_back<Seq,S>::type,NIL,NIL> VSub;
		mpl::for_each<T>( VSub() );
	}
};

template<typename Chooser, typename Seq> class SequenceRunner<Chooser, Seq, NIL, NIL>
{
public:
	template<typename S> void operator()(S)
	{
		typedef typename mpl::push_back<Seq,S>::type finalSeq;

		finalSeq f;

		Chooser c;
		c( f );
	}
};

template<typename Chooser, typename T, typename... Ts> class SequenceRunnerFrontend
{
public:
	typedef mpl::vector<> Seq;

	SequenceRunnerFrontend()
	{
		init();
	}


	void run( size_t id )
	{
		if( Chooser::id != id )
		{
			set( id );
		}
		Chooser::run( Chooser::object );
	}

	void set( size_t id )
	{
		Chooser::id = id;
		exec();
	}

private:
	void init()
	{
		Chooser::ids.clear();
		Chooser::init = true;
		exec();
		Chooser::init = false;
	}

	void exec()
	{
		typedef typename mpl::if_< mpl::bool_<(sizeof...(Ts) == 0)>, SequenceRunner<Chooser, Seq,NIL,NIL>, SequenceRunner<Chooser, Seq,Ts...,NIL> >::type VSub;
		mpl::for_each<T>( VSub() );
	}
};

#endif /* SEQUENCERRUNNER_H_ */
