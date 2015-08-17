//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "sequence/rtseq.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <string>

typedef RuntimeSeq<64> RtSeq;

typedef unsigned long long ull;

BOOST_AUTO_TEST_CASE( TestRtSeqSelector ) {
	BOOST_CHECK_EQUAL('G', nucl(RtSeq(10, "ACGTACGTAC")[2]));
	BOOST_CHECK_EQUAL('G', nucl(RtSeq(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[2]));
	BOOST_CHECK_EQUAL('G', nucl(RtSeq(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[16]));
	BOOST_CHECK_EQUAL('T', nucl(RtSeq(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[17]));
	BOOST_CHECK_EQUAL('A', nucl(RtSeq(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[18]));
	BOOST_CHECK_EQUAL('C', nucl(RtSeq(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[19]));

    BOOST_CHECK_EQUAL('C', nucl(RtSeq(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACCT")[15]));
    BOOST_CHECK_EQUAL('G', nucl(RtSeq(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACCT")[16]));
    BOOST_CHECK_EQUAL('C', nucl(RtSeq(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACCT")[31]));
    BOOST_CHECK_EQUAL('G', nucl(RtSeq(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACCT")[32]));
    BOOST_CHECK_EQUAL('T', nucl(RtSeq(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACCT")[63]));

    BOOST_CHECK_EQUAL('C', nucl(RtSeq(32, "ACGTACGTACACGTACGTACACGTACGTACAC")[15]));
    BOOST_CHECK_EQUAL('G', nucl(RtSeq(32, "ACGTACGTACACGTACGTACACGTACGTACAC")[16]));
    BOOST_CHECK_EQUAL('C', nucl(RtSeq(32, "ACGTACGTACACGTACGTACACGTACGTACAC")[31]));

    BOOST_CHECK_EQUAL('C', nucl(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACC")[15]));
    BOOST_CHECK_EQUAL('G', nucl(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACC")[16]));
    BOOST_CHECK_EQUAL('C', nucl(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACC")[32]));

}

BOOST_AUTO_TEST_CASE( TestRtSeqShiftLeft ) {
    BOOST_CHECK_EQUAL(RtSeq(5, "ACACA"), (RtSeq(5, "CACAC") << dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(5, "ACACC"), (RtSeq(5, "CACAC") << dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(5, "ACACG"), (RtSeq(5, "CACAC") << dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(5, "ACACT"), (RtSeq(5, "CACAC") << dignucl('T')));

    RtSeq s(10, "ACGTACGTAC");
	BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACA"), (s << dignucl('A')));
	BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACC"), (s << dignucl('C')));
	BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACG"), (s << dignucl('G')));
	BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACT"), (s << dignucl('T')));

	RtSeq s2(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
	BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACA"), (s2 << dignucl('A')));
	BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACC"), (s2 << dignucl('C')));
	BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACG"), (s2 << dignucl('G')));
	BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACT"), (s2 << dignucl('T')));

    RtSeq s2b(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGA"), (s2b << dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGC"), (s2b << dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGG"), (s2b << dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGT"), (s2b << dignucl('T')));

    RtSeq s3(32, "ACGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACA"), (s3 << dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACC"), (s3 << dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACG"), (s3 << dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACT"), s3 << dignucl('T'));

    RtSeq s4(33, "TACGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACA"), (s4 << dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACC"), (s4 << dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACG"), (s4 << dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACT"), s4 << dignucl('T'));
}


BOOST_AUTO_TEST_CASE( TestRtSeqShiftLeftThis ) {
    RtSeq s0(5, "CACAC");
    s0 <<= 'A';
    BOOST_CHECK_EQUAL(RtSeq(5, "ACACA"), s0);
    s0 <<= 'C';
    BOOST_CHECK_EQUAL(RtSeq(5, "CACAC"), s0);
    s0 <<= 'G';
    BOOST_CHECK_EQUAL(RtSeq(5, "ACACG"), s0);
    s0 <<= 'T';
    BOOST_CHECK_EQUAL(RtSeq(5, "CACGT"), s0);

    size_t l;
    RtSeq s(10, "ACGTACGTAC");
    l = 10;
    s <<= 'A';
    BOOST_CHECK_EQUAL(RtSeq(l, "CGTACGTACA"), s);
    s <<= 'C';
    BOOST_CHECK_EQUAL(RtSeq(l, "GTACGTACAC"), s);
    s <<= 'G';
    BOOST_CHECK_EQUAL(RtSeq(l, "TACGTACACG"), s);
    s <<= 'T';
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACACGT"), s);

    RtSeq s2(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    l = 60;
    s2 <<= 'A';
    BOOST_CHECK_EQUAL(RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACA"), s2);
    s2 <<= 'C';
    BOOST_CHECK_EQUAL(RtSeq(l, "GTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACAC"), s2);
    s2 <<= 'G';
    BOOST_CHECK_EQUAL(RtSeq(l, "TACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACG"), s2);
    s2 <<= 'T';
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGT"), s2);

    RtSeq s2b(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    l = 64;
    s2b <<= 'A';
    BOOST_CHECK_EQUAL(RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGA"), s2b);
    s2b <<= 'C';
    BOOST_CHECK_EQUAL(RtSeq(l, "GTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGAC"), s2b);
    s2b <<= 'G';
    BOOST_CHECK_EQUAL(RtSeq(l, "TACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGACG"), s2b);
    s2b <<= 'T';
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGACGT"), s2b);

    RtSeq s3(32, "ACGTACGTACACGTACGTACACGTACGTACAC");
    l = 32;
    s3 <<= 'C';
    BOOST_CHECK_EQUAL(RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACC"), s3);
    s3 <<= 'T';
    BOOST_CHECK_EQUAL(RtSeq(l, "GTACGTACACGTACGTACACGTACGTACACCT"), s3);
    s3 <<= 'A';
    BOOST_CHECK_EQUAL(RtSeq(l, "TACGTACACGTACGTACACGTACGTACACCTA"), s3);
    s3 <<= 'G';
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACACGTACGTACACGTACGTACACCTAG"), s3);

    RtSeq s4(33, "TACGTACGTACACGTACGTACACGTACGTACAC");
    l = 33;
    s4 <<= 'C';
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACGTACACGTACGTACACGTACGTACACC"), s4);
    s4 <<= 'T';
    BOOST_CHECK_EQUAL(RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACCT"), s4);
    s4 <<= 'A';
    BOOST_CHECK_EQUAL(RtSeq(l, "GTACGTACACGTACGTACACGTACGTACACCTA"), s4);
    s4 <<= 'G';
    BOOST_CHECK_EQUAL(RtSeq(l, "TACGTACACGTACGTACACGTACGTACACCTAG"), s4);
}


BOOST_AUTO_TEST_CASE( TestRtSeqShiftRight ) {
    BOOST_CHECK_EQUAL(RtSeq(5, "ACACA"), (RtSeq(5, "CACAC") >> dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(5, "CCACA"), (RtSeq(5, "CACAC") >> dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(5, "GCACA"), (RtSeq(5, "CACAC") >> dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(5, "TCACA"), (RtSeq(5, "CACAC") >> dignucl('T')));

	RtSeq s(10, "ACGTACGTAC");
	BOOST_CHECK_EQUAL(RtSeq(10, "AACGTACGTA"), (s >> dignucl('A')));
	BOOST_CHECK_EQUAL(RtSeq(10, "CACGTACGTA"), (s >> dignucl('C')));
	BOOST_CHECK_EQUAL(RtSeq(10, "GACGTACGTA"), (s >> dignucl('G')));
	BOOST_CHECK_EQUAL(RtSeq(10, "TACGTACGTA"), (s >> dignucl('T')));

	RtSeq s2(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
	BOOST_CHECK_EQUAL(RtSeq(60, "AACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA"), (s2 >> dignucl('A')));
	BOOST_CHECK_EQUAL(RtSeq(60, "CACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA"), (s2 >> dignucl('C')));
	BOOST_CHECK_EQUAL(RtSeq(60, "GACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA"), (s2 >> dignucl('G')));
	BOOST_CHECK_EQUAL(RtSeq(60, "TACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA"), (s2 >> dignucl('T')));

    RtSeq s2b(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    BOOST_CHECK_EQUAL(RtSeq(64, "AACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGC"), (s2b >> dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(64, "CACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGC"), (s2b >> dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(64, "GACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGC"), (s2b >> dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(64, "TACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGC"), (s2b >> dignucl('T')));

    RtSeq s3(32, "ACGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(32, "AACGTACGTACACGTACGTACACGTACGTACA"), (s3 >> dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(32, "CACGTACGTACACGTACGTACACGTACGTACA"), (s3 >> dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(32, "GACGTACGTACACGTACGTACACGTACGTACA"), (s3 >> dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(32, "TACGTACGTACACGTACGTACACGTACGTACA"), (s3 >> dignucl('T')));

    RtSeq s4(33, "ACGTACGTACACGTACGTACACGTACGTACACT");
    BOOST_CHECK_EQUAL(RtSeq(33, "AACGTACGTACACGTACGTACACGTACGTACAC"), (s4 >> dignucl('A')));
    BOOST_CHECK_EQUAL(RtSeq(33, "CACGTACGTACACGTACGTACACGTACGTACAC"), (s4 >> dignucl('C')));
    BOOST_CHECK_EQUAL(RtSeq(33, "GACGTACGTACACGTACGTACACGTACGTACAC"), (s4 >> dignucl('G')));
    BOOST_CHECK_EQUAL(RtSeq(33, "TACGTACGTACACGTACGTACACGTACGTACAC"), (s4 >> dignucl('T')));
}


BOOST_AUTO_TEST_CASE( TestRtSeqShiftRightThis ) {
    size_t l;
    RtSeq s(5, "CACAC");
    l = 5;
    s >>= ('A');
    BOOST_CHECK_EQUAL(RtSeq(l, "ACACA"), s);
    s >>= ('C');
    BOOST_CHECK_EQUAL(RtSeq(l, "CACAC"), s);
    s >>= ('G');
    BOOST_CHECK_EQUAL(RtSeq(l, "GCACA"), s);
    s >>= ('T');
    BOOST_CHECK_EQUAL(RtSeq(l, "TGCAC"), s);


    l = 9;
    s = RtSeq(l, "CGTACGTAC");
    s >>= ('A');
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACGTA"), s);
    s >>= ('C');
    BOOST_CHECK_EQUAL(RtSeq(l, "CACGTACGT"), s);
    s >>= ('G');
    BOOST_CHECK_EQUAL(RtSeq(l, "GCACGTACG"), s);
    s >>= ('T');
    BOOST_CHECK_EQUAL(RtSeq(l, "TGCACGTAC"), s);

    l = 59;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    s >>= ('A');
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA"), s);
    s >>= ('C');
    BOOST_CHECK_EQUAL(RtSeq(l, "CACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGT"), s);
    s >>= ('G');
    BOOST_CHECK_EQUAL(RtSeq(l, "GCACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACG"), s);
    s >>= ('T');
    BOOST_CHECK_EQUAL(RtSeq(l, "TGCACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTAC"), s);

    l = 63;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s >>= ('A');
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGC"), s);
    s >>= ('C');
    BOOST_CHECK_EQUAL(RtSeq(l, "CACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCG"), s);
    s >>= ('G');
    BOOST_CHECK_EQUAL(RtSeq(l, "GCACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACC"), s);
    s >>= ('T');
    BOOST_CHECK_EQUAL(RtSeq(l, "TGCACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), s);


    l = 31;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACAC");
    s >>= ('A');
    BOOST_CHECK_EQUAL(RtSeq(l, "ACGTACGTACACGTACGTACACGTACGTACA"), s);
    s >>= ('C');
    BOOST_CHECK_EQUAL(RtSeq(l, "CACGTACGTACACGTACGTACACGTACGTAC"), s);
    s >>= ('G');
    BOOST_CHECK_EQUAL(RtSeq(l, "GCACGTACGTACACGTACGTACACGTACGTA"), s);
    s >>= ('T');
    BOOST_CHECK_EQUAL(RtSeq(l, "TGCACGTACGTACACGTACGTACACGTACGT"), s);

    l = 32;
    s = RtSeq(l, "ACGTACGTACACGTACGTACACGTACGTACAC");
    s >>= ('A');
    BOOST_CHECK_EQUAL(RtSeq(l, "AACGTACGTACACGTACGTACACGTACGTACA"), s);
    s >>= ('C');
    BOOST_CHECK_EQUAL(RtSeq(l, "CAACGTACGTACACGTACGTACACGTACGTAC"), s);
    s >>= ('G');
    BOOST_CHECK_EQUAL(RtSeq(l, "GCAACGTACGTACACGTACGTACACGTACGTA"), s);
    s >>= ('T');
    BOOST_CHECK_EQUAL(RtSeq(l, "TGCAACGTACGTACACGTACGTACACGTACGT"), s);

    l = 33;
    s = RtSeq(l, "CACGTACGTACACGTACGTACACGTACGTACAC");
    s >>= ('A');
    BOOST_CHECK_EQUAL(RtSeq(l, "ACACGTACGTACACGTACGTACACGTACGTACA"), s);
    s >>= ('C');
    BOOST_CHECK_EQUAL(RtSeq(l, "CACACGTACGTACACGTACGTACACGTACGTAC"), s);
    s >>= ('G');
    BOOST_CHECK_EQUAL(RtSeq(l, "GCACACGTACGTACACGTACGTACACGTACGTA"), s);
    s >>= ('T');
    BOOST_CHECK_EQUAL(RtSeq(l, "TGCACACGTACGTACACGTACGTACACGTACGT"), s);
}



BOOST_AUTO_TEST_CASE( TestRtSeqStr ) {
    RtSeq s(10, "ACGTACGTAC");
	BOOST_CHECK_EQUAL("ACGTACGTAC", s.str());
	RtSeq s2(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
	BOOST_CHECK_EQUAL("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", s2.str());
    RtSeq s2b(64, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG", s2b.str());
    RtSeq s3(32, "TCGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL("TCGTACGTACACGTACGTACACGTACGTACAC", s3.str());
    RtSeq s4(33, "ACGTACGTACACGTACGTACACGTACGTACACT");
    BOOST_CHECK_EQUAL("ACGTACGTACACGTACGTACACGTACGTACACT", s4.str());
}

BOOST_AUTO_TEST_CASE( TestRtSeqHeadAndTail ) {
    size_t l;

    l = 5;
    RtSeq s(l, "GCATC");
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "CATC"), RtSeq(l - 1, s, 1)); // tail
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "GCAT"), RtSeq(l - 1, s)); // head

    l = 10;
    RtSeq s1(l, "CCGTACGTAC");
	BOOST_CHECK_EQUAL(RtSeq(l - 1, "CGTACGTAC"), RtSeq(l - 1, s1, 1)); // tail
	BOOST_CHECK_EQUAL(RtSeq(l - 1, "CCGTACGTA"), RtSeq(l - 1, s1)); // head

	l = 60;
    RtSeq s2(l, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), RtSeq(l - 1, s2, 1)); // tail
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA"), RtSeq(l - 1, s2)); // head

    l = 64;
    RtSeq s2b(l, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), RtSeq(l - 1, s2b, 1)); // tail
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGC"), RtSeq(l - 1, s2b)); // head

    l = 32;
    RtSeq s3(l, "TCGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "CGTACGTACACGTACGTACACGTACGTACAC"), RtSeq(l - 1, s3, 1)); // tail
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "TCGTACGTACACGTACGTACACGTACGTACA"), RtSeq(l - 1, s3)); // head

    l = 33;
    RtSeq s4(l, "GCGTACGTACACGTACGTACACGTACGTACACT");
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "CGTACGTACACGTACGTACACGTACGTACACT"), RtSeq(l - 1, s4, 1)); // tail
    BOOST_CHECK_EQUAL(RtSeq(l - 1, "GCGTACGTACACGTACGTACACGTACGTACAC"), RtSeq(l - 1, s4)); // head
}

BOOST_AUTO_TEST_CASE( TestRtSeqFromBiggerRtSeq ) {
    RtSeq s(64, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
	BOOST_CHECK_EQUAL("CCGTA", RtSeq(5, s).str());
	BOOST_CHECK_EQUAL("CCGTACGTAC", RtSeq(10, s).str());
	BOOST_CHECK_EQUAL("CCGTACGTACACGTAC", RtSeq(16, s).str());
	BOOST_CHECK_EQUAL("CCGTACGTACACGTACG", RtSeq(17, s).str());
	BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACAC", RtSeq(32, s).str());
	BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACACG", RtSeq(33, s).str());
	BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", RtSeq(60, s).str());
	BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG", RtSeq(64, s).str());
}

BOOST_AUTO_TEST_CASE( TestRtSeqFromSeq ) {
    Seq<5> s("CGTAC");
    Seq<10> s1("ACGTACGTAC");
    Seq<16> s2("CCCCGTACGTACGTAC");
    Seq<32> s3("ACGTACGTACACGTACGTACACGTACGTACAC");
    Seq<33> s4("GACGTACGTACACGTACGTACACGTACGTACAC");
    Seq<60> s5("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    Seq<64> s6("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");

    BOOST_CHECK_EQUAL(RtSeq(s, true).str(), s.str());
    BOOST_CHECK_EQUAL(RtSeq(s1, true).str(), s1.str());
    BOOST_CHECK_EQUAL(RtSeq(s2, true).str(), s2.str());
    BOOST_CHECK_EQUAL(RtSeq(s3, true).str(), s3.str());
    BOOST_CHECK_EQUAL(RtSeq(s4, true).str(), s4.str());
    BOOST_CHECK_EQUAL(RtSeq(s5, true).str(), s5.str());
    BOOST_CHECK_EQUAL(RtSeq(s6, true).str(), s6.str());
}

BOOST_AUTO_TEST_CASE( TestRtSeqToSeq ) {
    Seq<5> s("CGTAC");
    Seq<10> s1("ACGTACGTAC");
    Seq<16> s2("CCCCGTACGTACGTAC");
    Seq<32> s3("ACGTACGTACACGTACGTACACGTACGTACAC");
    Seq<33> s4("GACGTACGTACACGTACGTACACGTACGTACAC");
    Seq<60> s5("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    Seq<64> s6("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");

    BOOST_CHECK_EQUAL(RtSeq(5, "CGTAC").get_seq<5>(), s);
    BOOST_CHECK_EQUAL(RtSeq(10, "ACGTACGTAC").get_seq<10>(), s1);
    BOOST_CHECK_EQUAL(RtSeq(16, "CCCCGTACGTACGTAC").get_seq<16>(), s2);
    BOOST_CHECK_EQUAL(RtSeq(32, "ACGTACGTACACGTACGTACACGTACGTACAC").get_seq<32>(), s3);
    BOOST_CHECK_EQUAL(RtSeq(33, "GACGTACGTACACGTACGTACACGTACGTACAC").get_seq<33>(), s4);
    BOOST_CHECK_EQUAL(RtSeq(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC").get_seq<60>(), s5);
    BOOST_CHECK_EQUAL(RtSeq(64, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG").get_seq<64>(), s6);

}

BOOST_AUTO_TEST_CASE( TestRtSeqFromType ) {
    Sequence s("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    BOOST_CHECK_EQUAL("CCGTA", RtSeq(5, s).str());
    BOOST_CHECK_EQUAL("CCGTACGTAC", RtSeq(10, s).str());
    BOOST_CHECK_EQUAL("CCGTACGTACACGTAC", RtSeq(16, s).str());
    BOOST_CHECK_EQUAL("CCGTACGTACACGTACG", RtSeq(17, s).str());
    BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACAC", RtSeq(32, s).str());
    BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACACG", RtSeq(33, s).str());
    BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", RtSeq(60, s).str());
    BOOST_CHECK_EQUAL("CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG", RtSeq(64, s).str());

    BOOST_CHECK_EQUAL("GTA", RtSeq(3, s, 2).str());
    BOOST_CHECK_EQUAL("GTACGTACAC", RtSeq(10, s, 2).str());
    BOOST_CHECK_EQUAL("GTACGTACACGTACGT", RtSeq(16, s, 2).str());
    BOOST_CHECK_EQUAL("GTACGTACACGTACGTA", RtSeq(17, s, 2).str());
    BOOST_CHECK_EQUAL("GTACGTACACGTACGTACACGTACGTACACGT", RtSeq(32, s, 2).str());
    BOOST_CHECK_EQUAL("GTACGTACACGTACGTACACGTACGTACACGTA", RtSeq(33, s, 2).str());
    BOOST_CHECK_EQUAL("GTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", RtSeq(58, s, 2).str());
    BOOST_CHECK_EQUAL("GTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG", RtSeq(62, s, 2).str());
}

BOOST_AUTO_TEST_CASE( TestRtSeqPushBack ) {
    BOOST_CHECK_EQUAL(RtSeq(6, "CACACA"), (RtSeq(5, "CACAC").pushBack(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(6, "CACACC"), (RtSeq(5, "CACAC").pushBack(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(6, "CACACG"), (RtSeq(5, "CACAC").pushBack(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(6, "CACACT"), (RtSeq(5, "CACAC").pushBack(dignucl('T'))));

    RtSeq s(9, "CGTACGTAC");
    BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACA"), (s.pushBack(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACC"), (s.pushBack(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACG"), (s.pushBack(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(10, "CGTACGTACT"), (s.pushBack(dignucl('T'))));

    RtSeq s2(59, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACA"), (s2.pushBack(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACC"), (s2.pushBack(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACG"), (s2.pushBack(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(60, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACT"), (s2.pushBack(dignucl('T'))));

    RtSeq s2b(63, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGA"), (s2b.pushBack(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGC"), (s2b.pushBack(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGG"), (s2b.pushBack(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(64, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGT"), (s2b.pushBack(dignucl('T'))));

    RtSeq s3(31, "CGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACA"), (s3.pushBack(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACC"), (s3.pushBack(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACG"), (s3.pushBack(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(32, "CGTACGTACACGTACGTACACGTACGTACACT"), (s3.pushBack(dignucl('T'))));

    RtSeq s4(32, "ACGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACA"), (s4.pushBack(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACC"), (s4.pushBack(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACG"), (s4.pushBack(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(33, "ACGTACGTACACGTACGTACACGTACGTACACT"), (s4.pushBack(dignucl('T'))));

    RtSeq s5(33, "CACGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(34, "CACGTACGTACACGTACGTACACGTACGTACACA"), (s5.pushBack(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(34, "CACGTACGTACACGTACGTACACGTACGTACACC"), (s5.pushBack(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(34, "CACGTACGTACACGTACGTACACGTACGTACACG"), (s5.pushBack(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(34, "CACGTACGTACACGTACGTACACGTACGTACACT"), (s5.pushBack(dignucl('T'))));
}


BOOST_AUTO_TEST_CASE( TestRtSeqPushBackThis ) {
    size_t l;
    RtSeq s(5, "CACAC");
    l = 5;
    s.pushBackThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CACACA"), s);
    s.pushBackThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CACACAC"), s);
    s.pushBackThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "CACACACG"), s);
    s.pushBackThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "CACACACGT"), s);


    l = 9;
    s = RtSeq(l, "CGTACGTAC");
    s.pushBackThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CGTACGTACA"), s);
    s.pushBackThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CGTACGTACAC"), s);
    s.pushBackThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "CGTACGTACACG"), s);
    s.pushBackThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "CGTACGTACACGT"), s);

    l = 59;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    s.pushBackThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACA"), s);
    s.pushBackThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushBackThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACG"), s);
    s.pushBackThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGT"), s);

    l = 63;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushBackThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGA"), s);
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushBackThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGC"), s);
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushBackThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGG"), s);
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushBackThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCGT"), s);


    l = 31;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACAC");
    s.pushBackThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CGTACGTACACGTACGTACACGTACGTACACA"), s);
    s.pushBackThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CGTACGTACACGTACGTACACGTACGTACACAC"), s);
    s.pushBackThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "CGTACGTACACGTACGTACACGTACGTACACACG"), s);
    s.pushBackThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "CGTACGTACACGTACGTACACGTACGTACACACGT"), s);

    l = 32;
    s = RtSeq(l, "ACGTACGTACACGTACGTACACGTACGTACAC");
    s.pushBackThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "ACGTACGTACACGTACGTACACGTACGTACACA"), s);
    s.pushBackThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "ACGTACGTACACGTACGTACACGTACGTACACAC"), s);
    s.pushBackThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "ACGTACGTACACGTACGTACACGTACGTACACACG"), s);
    s.pushBackThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "ACGTACGTACACGTACGTACACGTACGTACACACGT"), s);

    l = 33;
    s = RtSeq(l, "CACGTACGTACACGTACGTACACGTACGTACAC");
    s.pushBackThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CACGTACGTACACGTACGTACACGTACGTACACA"), s);
    s.pushBackThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CACGTACGTACACGTACGTACACGTACGTACACAC"), s);
    s.pushBackThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "CACGTACGTACACGTACGTACACGTACGTACACACG"), s);
    s.pushBackThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "CACGTACGTACACGTACGTACACGTACGTACACACGT"), s);
}


BOOST_AUTO_TEST_CASE( TestRtSeqPushFront ) {
    BOOST_CHECK_EQUAL(RtSeq(6, "ACACAC"), (RtSeq(5, "CACAC").pushFront(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(6, "CCACAC"), (RtSeq(5, "CACAC").pushFront(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(6, "GCACAC"), (RtSeq(5, "CACAC").pushFront(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(6, "TCACAC"), (RtSeq(5, "CACAC").pushFront(dignucl('T'))));

    RtSeq s(9, "CGTACGTAC");
    BOOST_CHECK_EQUAL(RtSeq(10, "ACGTACGTAC"), (s.pushFront(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(10, "CCGTACGTAC"), (s.pushFront(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(10, "GCGTACGTAC"), (s.pushFront(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(10, "TCGTACGTAC"), (s.pushFront(dignucl('T'))));

    RtSeq s2(59, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL(RtSeq(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), (s2.pushFront(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(60, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), (s2.pushFront(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(60, "GCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), (s2.pushFront(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(60, "TCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), (s2.pushFront(dignucl('T'))));

    RtSeq s2b(63, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    BOOST_CHECK_EQUAL(RtSeq(64, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), (s2b.pushFront(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(64, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), (s2b.pushFront(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(64, "GCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), (s2b.pushFront(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(64, "TCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), (s2b.pushFront(dignucl('T'))));

    RtSeq s3(31, "CGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(32, "ACGTACGTACACGTACGTACACGTACGTACAC"), (s3.pushFront(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(32, "CCGTACGTACACGTACGTACACGTACGTACAC"), (s3.pushFront(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(32, "GCGTACGTACACGTACGTACACGTACGTACAC"), (s3.pushFront(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(32, "TCGTACGTACACGTACGTACACGTACGTACAC"), (s3.pushFront(dignucl('T'))));

    RtSeq s4(32, "ACGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(33, "AACGTACGTACACGTACGTACACGTACGTACAC"), (s4.pushFront(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(33, "CACGTACGTACACGTACGTACACGTACGTACAC"), (s4.pushFront(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(33, "GACGTACGTACACGTACGTACACGTACGTACAC"), (s4.pushFront(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(33, "TACGTACGTACACGTACGTACACGTACGTACAC"), (s4.pushFront(dignucl('T'))));

    RtSeq s5(33, "CACGTACGTACACGTACGTACACGTACGTACAC");
    BOOST_CHECK_EQUAL(RtSeq(34, "ACACGTACGTACACGTACGTACACGTACGTACAC"), (s5.pushFront(dignucl('A'))));
    BOOST_CHECK_EQUAL(RtSeq(34, "CCACGTACGTACACGTACGTACACGTACGTACAC"), (s5.pushFront(dignucl('C'))));
    BOOST_CHECK_EQUAL(RtSeq(34, "GCACGTACGTACACGTACGTACACGTACGTACAC"), (s5.pushFront(dignucl('G'))));
    BOOST_CHECK_EQUAL(RtSeq(34, "TCACGTACGTACACGTACGTACACGTACGTACAC"), (s5.pushFront(dignucl('T'))));
}


BOOST_AUTO_TEST_CASE( TestRtSeqPushFrontThis ) {
    size_t l;
    RtSeq s(5, "CACAC");
    l = 5;
    s.pushFrontThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "ACACAC"), s);
    s.pushFrontThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CACACAC"), s);
    s.pushFrontThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "GCACACAC"), s);
    s.pushFrontThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "TGCACACAC"), s);


    l = 9;
    s = RtSeq(l, "CGTACGTAC");
    s.pushFrontThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "ACGTACGTAC"), s);
    s.pushFrontThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CACGTACGTAC"), s);
    s.pushFrontThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "GCACGTACGTAC"), s);
    s.pushFrontThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "TGCACGTACGTAC"), s);

    l = 59;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    s.pushFrontThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), s);
    s.pushFrontThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), s);
    s.pushFrontThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "GCACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), s);
    s.pushFrontThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "TGCACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"), s);

    l = 63;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushFrontThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), s);
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushFrontThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "CCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), s);
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushFrontThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "GCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), s);
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG");
    s.pushFrontThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "TCGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACCGCG"), s);


    l = 31;
    s = RtSeq(l, "CGTACGTACACGTACGTACACGTACGTACAC");
    s.pushFrontThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "ACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "GCACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "TGCACGTACGTACACGTACGTACACGTACGTACAC"), s);

    l = 32;
    s = RtSeq(l, "ACGTACGTACACGTACGTACACGTACGTACAC");
    s.pushFrontThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "AACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CAACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "GCAACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "TGCAACGTACGTACACGTACGTACACGTACGTACAC"), s);

    l = 33;
    s = RtSeq(l, "CACGTACGTACACGTACGTACACGTACGTACAC");
    s.pushFrontThis('A');
    BOOST_CHECK_EQUAL(RtSeq(l + 1, "ACACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('C');
    BOOST_CHECK_EQUAL(RtSeq(l + 2, "CACACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('G');
    BOOST_CHECK_EQUAL(RtSeq(l + 3, "GCACACGTACGTACACGTACGTACACGTACGTACAC"), s);
    s.pushFrontThis('T');
    BOOST_CHECK_EQUAL(RtSeq(l + 4, "TGCACACGTACGTACACGTACGTACACGTACGTACAC"), s);
}

BOOST_AUTO_TEST_CASE( TestRtSeqNull ) {
    RtSeq s(0,"");
	BOOST_CHECK_EQUAL("", s.str());
}


BOOST_AUTO_TEST_CASE( TestRtSeqAddSymbolForNullValue ) {
    RtSeq s1(1, "G");
    RtSeq s2 = (s1 << 'A');
    RtSeq s3(1, "A");
	BOOST_CHECK_EQUAL(s3, s2);
}


BOOST_AUTO_TEST_CASE( TestRtSeqComplex ) {
	Sequence s1("ACAAA");
	Sequence s2("CAAAC");
	BOOST_CHECK_EQUAL((!(RtSeq(4, !s1))).str(), RtSeq(4, s2).str());
	BOOST_CHECK_EQUAL(!(RtSeq(4, !s1)), RtSeq(4, s2));
}

BOOST_AUTO_TEST_CASE( TestRtSeqFromCharArray ) {
	std::string s = "ACGTACGTAC";
	BOOST_CHECK_EQUAL("ACGTACGTAC", RtSeq(10, s.c_str()).str());
}

BOOST_AUTO_TEST_CASE( TestRtSeqReverseComplement ) {
    RtSeq s(10, "ACGTACGTAC");
    BOOST_CHECK_EQUAL("GTACGTACGT", (!s).str());
    RtSeq s1(9, "CGTACGTAC");
    BOOST_CHECK_EQUAL("GTACGTACG", (!s1).str());
    RtSeq s2(60, "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL("GTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGT", (!s2).str());
    RtSeq s2b(64, "TGCAACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL("GTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTTGCA", (!s2b).str());
    RtSeq s3(32, "CGTACGTACACGTACGTACACGTACGTACACG");
    BOOST_CHECK_EQUAL("CGTGTACGTACGTGTACGTACGTGTACGTACG", (!s3).str());
    RtSeq s4(33, "ACGTACGTACACGTACGTACACGTACGTACACG");
    BOOST_CHECK_EQUAL("CGTGTACGTACGTGTACGTACGTGTACGTACGT", (!s4).str());
}

BOOST_AUTO_TEST_CASE( TestRtSeq16 ) {
    RtSeq s(16, "AAAAAAAAAAAAAAAA");
	BOOST_CHECK_EQUAL(s << 'C', RtSeq(16, "AAAAAAAAAAAAAAAC"));
}

BOOST_AUTO_TEST_CASE( TestRtSeq16_2 ) {
    RtSeq s(16, "TTTTTTTTTTTTTTTT");
	BOOST_CHECK_EQUAL(RtSeq(16, "TTTTTTTTTTTTTTTA"), s << 'A');
}

BOOST_AUTO_TEST_CASE( TestRtSeqFirstLast ) {
    RtSeq s1(7, "ACGTACT");
	BOOST_CHECK_EQUAL(0, s1.first());
	BOOST_CHECK_EQUAL(3, s1.last());
	RtSeq s2(7, "TTTTTTT");
	BOOST_CHECK_EQUAL(3, s2.first());
	BOOST_CHECK_EQUAL(3, s2.last());
}
