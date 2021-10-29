/**
 * (c) Camille Scott, 2019
 * File   : processors.hh
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 15.10.2019
 *
 * Sequence-parsing driver base classes. These classes manage a 
 * sequence parser and delegate processing to a process_sequence method.
 * The InserterProcessor is a generic processor that parses reads
 * and repeatedly calls insert_sequence on any class that has the method.
 *
 */

#ifndef GOETIA_PROCESSORS_HH
#define GOETIA_PROCESSORS_HH

#include <tuple>
#include <memory>
#include <fstream>
#include <iostream>
#include <string>

#include "goetia/goetia.hh"
#include "goetia/metrics.hh"
#include "goetia/parsing/parsing.hh"
#include "goetia/parsing/readers.hh"
#include "goetia/sequences/exceptions.hh"


namespace goetia {


/**
 * @Synopsis  CRTP base class for generic sequence processing. Uses
 *            the given sequence parsing type to parse sequences and
 *            pass them to a process_sequence method that is implemented
 *            by the derived class. Sequences are processed in chunks
 *            corresponding to interval: the interval timing is incremented
 *            based on the return value of process_sequence.
 *
 * @tparam Derived    CRTP derived class.
 * @tparam ParserType Sequencing parsing type.
 */
template <class Derived,
          class ParserType = parsing::FastxParser<>>
class FileProcessor {

protected:

    metrics::IntervalCounter timer;
    metrics::Gauge           _n_sequences;
    bool                     _verbose;

public:

    typedef typename ParserType::alphabet alphabet;

    FileProcessor(uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL,
                  bool     verbose  = false)
        : timer(interval),
          _n_sequences{"timing", "n_sequences"},
          _verbose(verbose)
    {
        
    }

    /**
     * @Synopsis  Process a split pair of paired-end sequence files until
     *            all sequences are consumed.
     *
     * @Param left_filename    Left/R1 filename.
     * @Param right_filename   Right/R2 filename.
     * @Param min_length       Filter sequences under this length;
     *                         if 0 (default), do no filter.
     * @Param force_name_match Force left and right sequence names to follow
     *                         standard naming conventions.
     *
     * @Returns                Number of sequences processed, time passed.
     */
    std::tuple<uint64_t, uint64_t> process(const std::string& left_filename,
                     const std::string& right_filename,
                     bool strict = false,
                     uint32_t min_length=0,
                     bool force_name_match=false) {
        auto reader = parsing::SplitPairedReader<ParserType>::build(left_filename,
                                                                    right_filename,
                                                                    strict,
                                                                    min_length,
                                                                    force_name_match);
        return process(reader);
    }

    std::tuple<uint64_t, uint64_t> process(std::shared_ptr<parsing::SplitPairedReader<ParserType>>& reader) {
        uint64_t n_sequences, time_total;
        bool remaining = true;
        while(1) {
            std::tie(n_sequences, time_total, remaining) = advance(reader);
        }

        return {n_sequences, time_total};
    }

    /**
     * @Synopsis  Process a single-ended sequence file until all sequences
     *            are consumed.
     *
     * @Param filename File to process.
     *
     * @Returns   Number of sequences processed, time passed.
     */
    std::tuple<uint64_t, uint64_t> process(std::string const &filename,
                     bool strict = false,
                     uint32_t min_length = 0) {
        auto reader  = ParserType::build(filename, strict, min_length);
        return process(reader);
    }

    std::tuple<uint64_t, uint64_t> process(std::shared_ptr<ParserType>& reader) {
        uint64_t n_sequences = 0, time_total = 0;
        bool remaining = true;
        while(remaining) {
            std::tie(n_sequences, time_total, remaining) = advance(reader);
        }

        return {n_sequences, time_total};
    }

    /**
     * @Synopsis     Default paired-end sequence processing implementation:
     *               just consume both.
     *
     * @Param bundle ReadBundle containing the two sequences.
     */
    uint64_t process_sequence(parsing::RecordPair& pair) {
        uint64_t time_passed = 0;

        if (pair.first) {
            time_passed += derived().process_sequence(pair.first.value());
        }
        if (pair.second) {
            time_passed += derived().process_sequence(pair.second.value());
        }

        return time_passed;
    }

    template<typename ReaderType>
    auto handle_next(ReaderType& reader)
    -> std::optional<typename ReaderType::value_type> {
        try {
            auto record = reader.next();
            return record;
        } catch (InvalidCharacterException &e) {
            if (_verbose) {
                std::cerr << "WARNING: Bad sequence encountered at "
                          << this->n_sequences()
                          << ", exception was "
                          << e.what() << std::endl;
            }
            return {};
        }  catch (parsing::InvalidRead& e) {
            if (_verbose) {
                std::cerr << "WARNING: Invalid sequence encountered at "
                          << this->n_sequences()
                          << ", exception was "
                          << e.what() << std::endl;
            }
            return {};
        }
    }

    /**
     * @Synopsis  Consume sequences for the next [INTERVAL].
     *
     * @Param reader The reader to consume from.
     *
     * @Returns   Tuple containing <total seqs processed, total time passed, whether sequences remain>
     */
    std::tuple<uint64_t, uint64_t, bool> advance(std::shared_ptr<parsing::SplitPairedReader<ParserType>>& reader) {
        std::optional<parsing::RecordPair> bundle;
        while (!reader->is_complete()) {
            bundle = handle_next(*reader);
            
            if (!bundle) {
                continue;
            }

            uint64_t time_passed = derived().process_sequence(bundle.value());
            int _bundle_count = (bool)bundle.value().first + (bool)bundle.value().second;
            _n_sequences += _bundle_count;

            if (timer.poll(time_passed)) {
                return {_n_sequences, timer.total(), true};
            }
        }

        return {_n_sequences, timer.total(), false};
    }

    /**
     * @Synopsis  Consume sequences for the next [INTERVAL].
     *
     * @Param parser The reader to consume from.
     *
     * @Returns   Tuple containing <total seqs processed, total time passed, whether sequences remain>
     */
    std::tuple<uint64_t, uint64_t, bool> advance(std::shared_ptr<ParserType>& parser) {

        std::optional<parsing::Record> record;
        // Iterate through the reads and consume their k-mers.
        while (!parser->is_complete()) {
            record = handle_next(*parser);
            
            if (!record) {
                continue;
            }

            uint64_t time_passed = derived().process_sequence(record.value());
            ++_n_sequences;

            if (timer.poll(time_passed)) {
                return {_n_sequences, timer.total(), true};
            }
            
        }

        return {_n_sequences, timer.total(), false};
    }

    /**
     * @Synopsis  Number of sequences / reads consumed.
     *
     * @Returns   
     */
    uint64_t n_sequences() const {
        return _n_sequences;
    }

    uint64_t time_elapsed() const {
        return timer.total();
    }

    uint64_t interval() const {
        return timer.interval;
    }

    template<typename... Args>
    static std::shared_ptr<Derived> build(Args&&... args,
                                          uint64_t interval,
                                          bool     verbose) {
        return Derived::build(std::forward<Args>(args)...,
                              interval,
                              verbose);
    }

private:

    friend Derived;

    Derived& derived() {
        return *static_cast<Derived*>(this);
    }

    const Derived& derived() const {
        return *static_cast<const Derived*>(this);
    }

};


/**
 * @Synopsis  Generic processor for passing reads to a class
 *            with an `insert_sequence` method.
 *
 * @tparam InserterType Class with insert_sequence.
 * @tparam ParserType   Sequence parser type.
 */
template <class InserterType,
          class ParserType = parsing::FastxParser<>>
class InserterProcessor : public FileProcessor<InserterProcessor<InserterType, ParserType>,
                                               ParserType> {

protected:

    std::shared_ptr<InserterType> inserter;

    typedef FileProcessor<InserterProcessor<InserterType, ParserType>,
                          ParserType> Base;

public:

    using Base::process_sequence;
    typedef typename Base::alphabet alphabet;
    
    InserterProcessor(std::shared_ptr<InserterType> inserter,
                      uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL,
                      bool     verbose  = false)
        : Base(interval, verbose),
          inserter(inserter)
    {
    }

    uint64_t process_sequence(const parsing::Record& sequence) {
        try {
            return inserter->insert_sequence(sequence.sequence);
        } catch (SequenceLengthException &e) {
            if (this->_verbose) {
                std::cerr << "WARNING: Skipped sequence that was too short: sequence number "
                          << this->n_sequences() << " with sequence "
                          << sequence.sequence 
                          << std::endl;
            }
            return 0;
        } catch (InvalidCharacterException& e) {
            return 0;
        } catch (std::exception &e) {
            std::cerr << "ERROR: Exception thrown at " << this->n_sequences()
                      << " with msg: " << e.what()
                      << ". Sequence was: " << sequence.sequence
                      <<  std::endl;
            throw;
        }
    }

    void report() {

    }

    
    static auto build(std::shared_ptr<InserterType> inserter,
               uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL)
    -> std::shared_ptr<InserterProcessor<InserterType, ParserType>> {
        return std::make_shared<InserterProcessor<InserterType, ParserType>>(inserter, interval);
    }

};


/**
 * @Synopsis  Generic processor for passing reads to a class
 *            with a `filter_sequence` method.
 *
 * @tparam FilterType Class with filter_sequence.
 * @tparam ParserType Sequence parser type.
 */
template <class FilterType,
          class ParserType = parsing::FastxParser<>>
class FilterProcessor : public FileProcessor<FilterProcessor<FilterType, ParserType>,
                                             ParserType> {

protected:

    std::shared_ptr<FilterType> filter;
    std::ofstream               _output_stream; 
    metrics::Gauge              _n_passed;

    typedef FileProcessor<FilterProcessor<FilterType, ParserType>,
                          ParserType> Base;

public:

    using Base::process_sequence;
    typedef typename Base::alphabet alphabet;
    
    FilterProcessor(std::shared_ptr<FilterType> filter,
                    const std::string           output_filename,
                    uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL,
                    bool     verbose  = false)
        : Base(interval, verbose),
          filter(filter),
          _output_stream(output_filename.c_str()),
          _n_passed{"timing", "n_passed"}
    {
    }

    ~FilterProcessor() {
        _output_stream.close();
    }

    uint64_t process_sequence(const parsing::Record& sequence) {
        bool passed = false;
        uint64_t time_taken = 0;
        try {
            std::tie(passed, time_taken) = filter->filter_sequence(sequence.sequence);
        } catch (SequenceLengthException &e) {
            if (this->_verbose) {
                std::cerr << "WARNING: Skipped sequence that was too short: read "
                          << this->n_sequences() << " with sequence "
                          << sequence.sequence 
                          << std::endl;
            }
            return 0;
        } catch (InvalidCharacterException& e) {
            return 0;
        } catch (std::exception &e) {
            std::cerr << "ERROR: Exception thrown at " << this->n_sequences()
                      << " with msg: " << e.what()
                      <<  std::endl;
            throw e;
        }

        if (passed) {
            ++_n_passed;
            sequence.write_fastx(_output_stream);
        }

        return time_taken;
    }

    void report() {

    }

    uint64_t n_passed() const {
        return _n_passed;
    }

    static auto build(std::shared_ptr<FilterType> filter,
                      const std::string&          output_filename,
                      uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL,
                      bool verbose      = false)
    -> std::shared_ptr<FilterProcessor<FilterType, ParserType>> {

        return std::make_shared<FilterProcessor<FilterType,
                                                ParserType>>(filter,
                                                             output_filename,
                                                             interval,
                                                             verbose);
    }

};


/**
 * @Synopsis  Merges split-paired reads to single stream.
 *
 * @tparam ParserType Sequence parser type.
 */
template <class ParserType = parsing::FastxParser<>>
class MergeProcessor : public FileProcessor<MergeProcessor<ParserType>,
                                             ParserType> {

protected:

    std::ofstream               _output_stream; 
    metrics::Gauge              _n_passed;

    typedef FileProcessor<MergeProcessor<ParserType>,
                          ParserType> Base;

public:

    using Base::process_sequence;
    typedef typename Base::alphabet alphabet;
    
    MergeProcessor(const std::string           output_filename,
                   uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL,
                   bool     verbose  = false)
        : Base(interval, verbose),
          _output_stream(output_filename.c_str()),
          _n_passed{"timing", "n_passed"}
    {
    }

    ~MergeProcessor() {
        _output_stream.close();
    }

    uint64_t process_sequence(const parsing::Record& sequence) {
        sequence.write_fastx(_output_stream);

        return sequence.sequence.size();
    }

    void report() {

    }

    uint64_t n_passed() const {
        return _n_passed;
    }

    static auto build(const std::string&          output_filename,
                      uint64_t interval = metrics::IntervalCounter::DEFAULT_INTERVAL,
                      bool verbose      = false)
    -> std::shared_ptr<MergeProcessor<ParserType>> {

        return std::make_shared<MergeProcessor<ParserType>>(output_filename,
                                                             interval,
                                                             verbose);
    }

};



} //namespace goetia
#endif
