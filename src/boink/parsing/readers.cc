/**
 * (c) Camille Scott, 2019
 * File   : readers.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 21.10.2019
 */

#include "boink/parsing/readers.hh"

#include <fstream>
#include <string>
#include <utility>

// ignore warnings from seqan
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"

#pragma GCC diagnostic pop


namespace boink
{

namespace parsing
{


template<typename SeqIO>
SequenceReader<SeqIO>::SequenceReader(std::unique_ptr<SeqIO> pf)
{
    _parser = std::move(pf);
}


template<typename SeqIO>
SequenceReader<SeqIO>::SequenceReader(SequenceReader<SeqIO>& other)
{
    _parser = std::move(other._parser);
}


template<typename SeqIO>
SequenceReader<SeqIO>&
SequenceReader<SeqIO>::operator=(SequenceReader<SeqIO>& other) {
    _parser = std::move(other._parser);
    return *this;
}


template<typename SeqIO>
SequenceReader<SeqIO>::SequenceReader(SequenceReader<SeqIO>&&) noexcept {}


template<typename SeqIO>
std::shared_ptr<SequenceReader<SeqIO>>
SequenceReader<SeqIO>::build(const std::string& filename) {
    return std::make_shared<SequenceReader<SeqIO>>(
               std::move(std::make_unique<SeqIO>(filename))
           );
}


template<typename SeqIO>
Record SequenceReader<SeqIO>::next()
{
    return _parser->next();
}


template<typename SeqIO>
size_t SequenceReader<SeqIO>::num_parsed() const
{
    return _parser->num_parsed();
}

template<typename SeqIO>
bool SequenceReader<SeqIO>::is_complete() const {
    return _parser->is_complete();
}


FastxParser::FastxParser(const std::string& infile)
    : _filename(infile),
      _spin_lock(0),
      _num_parsed(0),
      _have_qualities(false),
      _is_complete(false)
{
    _fp = gzopen(_filename.c_str(), "r");
    _kseq = kseq_init(_fp);

    __asm__ __volatile__ ("" ::: "memory");
}

FastxParser::FastxParser()
    : FastxParser("-")
{
}


FastxParser::FastxParser(FastxParser& other) 
    : FastxParser(other._filename)
{   
}


FastxParser::FastxParser(FastxParser&& other) noexcept 
    : _filename(std::move(other._filename)),
      _spin_lock(other._spin_lock),
      _num_parsed(other._num_parsed),
      _have_qualities(other._have_qualities),
      _fp(other._fp),
      _kseq(other._kseq),
      _is_complete(other._is_complete)
{
    other._is_complete = true;
}


FastxParser::~FastxParser()
{
    kseq_destroy(_kseq);
    gzclose(_fp);
}


size_t FastxParser::num_parsed() const
{
    return _num_parsed;
}

bool FastxParser::is_complete() const {
    return _is_complete;
}


Record FastxParser::next()
{
    Record record;
    
    while (!__sync_bool_compare_and_swap(&_spin_lock, 0, 1));

    int stat = kseq_read(_kseq);
    if (stat >= 0) {
        record.sequence.assign(_kseq->seq.s, _kseq->seq.l);
        record.name.assign(_kseq->name.s, _kseq->name.l);
        if (_kseq->qual.l) {
            record.quality.assign(_kseq->qual.s, _kseq->qual.l);
            if (_num_parsed == 0) {
                _have_qualities = true;
            }
        }
        _num_parsed++;
    }

    __asm__ __volatile__ ("" ::: "memory");
    _spin_lock = 0;

    if (stat == -1) {
        _is_complete = true;
        throw NoMoreReadsAvailable();
    }

    if (stat == -2) {
        throw InvalidRead("Sequence and quality lengths differ");
    }

    if (stat == -3) {
        throw BoinkFileException("Error reading stream.");
    }

    return record;
}


// All template instantiations used in the codebase must be declared here.


}

}

template class boink::parsing::SequenceReader<boink::parsing::FastxParser>;
template class boink::parsing::SplitPairedReader<boink::parsing::FastxParser>;

// vim: set ft=cpp sts=4 sw=4 tw=80:

