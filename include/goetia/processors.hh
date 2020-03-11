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
#include "goetia/parsing/parsing.hh"
#include "goetia/parsing/readers.hh"
#include "goetia/events.hh"
#include "goetia/event_types.hh"

#include "goetia/sequences/exceptions.hh"


namespace goetia {

struct DEFAULT_INTERVALS {

    constexpr static unsigned long FINE = 10000;
    constexpr static unsigned long MEDIUM = 100000;
    constexpr static unsigned long COARSE = 1000000;

};



/**
 * @Synopsis  Bounded counter. Returns True when interval is
 *            reached and resets to zero.
 */
class IntervalCounter {

public:

    uint64_t interval;
    uint64_t counter;

    IntervalCounter(uint64_t interval)
        : interval(interval),
          counter(0)
    {
    }

    bool poll(uint64_t incr=1) {
        counter += incr;
        if (counter >= interval) {
            counter = 0;
            return true;
        } else {
            return false;
        }
    }

    friend inline std::ostream& operator<< (std::ostream& os, const IntervalCounter& counter) {
       os  << "IntervalCounter<interval=" << counter.interval 
           << ", counter=" << counter.counter << ">";
       return os;
    }
};


/**
 * @Synopsis  Reports the current sequence-processing interval.
 *            The intervals themselves are defined by the user
 *            and given to the processor; the processor notifies
 *            all its registered listeners with an interval_state
 *            stored on a TimeIntervalEvent when an interval
 *            is ticked.
 */
struct interval_state {
	bool fine;
	bool medium;
	bool coarse;
	bool end;

	interval_state()
		: fine(false), medium(false), coarse(false), end(false)
	{
	}

	interval_state(bool fine, bool medium, bool coarse, bool end)
		: fine(fine), medium(medium), coarse(coarse), end(end)
	{
	}
};


/**
 * @Synopsis  CRTP base class for generic sequence processing. Uses
 *            the given sequence parsing type to parse sequences and
 *            pass them to a process_sequence method that is implemented
 *            by the derived class. Notifies registered listeners when
 *            the given FINE, MEDIUM, and COARSE intervals of number
 *            of parsed sequences are reached.
 *
 * @tparam Derived    CRTP derived class.
 * @tparam ParserType Sequencing parsing type.
 */
template <class Derived,
          class ParserType = parsing::FastxParser<>>
class FileProcessor : public events::EventNotifier {

protected:

    std::array<IntervalCounter, 3> counters;
    uint64_t                       _n_reads;
    uint64_t                       _n_skipped;
    bool                           _verbose;

    bool _ticked(interval_state tick) {
        return tick.fine || tick.medium || tick.coarse || tick.end;
    }

    /**
     * @Synopsis  Increment all interval counters by n_ticks and notify
     *            listeners of interval reached.
     *
     * @Param     n_ticks Number of sequences to increment.
     * 
     * @Returns   An interval_state reporting whether an interval is reached.
     */
    interval_state _notify_tick(uint64_t n_ticks) {
        interval_state result;

        if (counters[0].poll(n_ticks)) {
             //std::cerr << "processed " << _n_reads << " sequences." << std::endl;               
             derived().report();
             auto event = std::make_shared<events::TimeIntervalEvent>();
             event->level = events::TimeIntervalEvent::FINE;
             event->t = _n_reads;
             notify(event);
             result.fine = true;
        }
        if (counters[1].poll(n_ticks)) {
             auto event = std::make_shared<events::TimeIntervalEvent>();
             event->level = events::TimeIntervalEvent::MEDIUM;
             event->t = _n_reads;
             notify(event);
             result.medium = true;
        }
        if (counters[2].poll(n_ticks)) {
             auto event = std::make_shared<events::TimeIntervalEvent>();
             event->level = events::TimeIntervalEvent::COARSE;
             event->t = _n_reads;
             notify(event);
             result.coarse = true;
        }

        result.end = false;
        return result;
    }

    /**
     * @Synopsis  Notifiy listeners that parsing is complete.
     */
    void _notify_stop() {
        auto event = std::make_shared<events::TimeIntervalEvent>();
        event->level = events::TimeIntervalEvent::END;
        event->t = _n_reads;
        notify(event);
    }


public:

    using events::EventNotifier::register_listener;
    using events::EventNotifier::notify;

    typedef typename ParserType::alphabet alphabet;

    FileProcessor(uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                  uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                  uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE,
                  bool     verbose         = false)
        :  events::EventNotifier(),
           counters { fine_interval, 
                      medium_interval,
                      coarse_interval },
          _n_reads(0),
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
     * @Returns                Number of sequences processed.
     */
    uint64_t process(const std::string& left_filename,
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

    uint64_t process(std::shared_ptr<parsing::SplitPairedReader<ParserType>>& reader) {
        while(1) {
            auto state = advance(reader);
            if (state.end) {
                break;
            }
        }

        return _n_reads;
    }

    /**
     * @Synopsis  Process a single-ended sequence file until all sequences
     *            are consumed.
     *
     * @Param filename File to process.
     *
     * @Returns   Number of sequences consumed.
     */
    uint64_t process(std::string const &filename,
                     bool strict = false,
                     uint32_t min_length = 0) {
        auto reader  = ParserType::build(filename, strict, min_length);
        return process(reader);
    }

    uint64_t process(std::shared_ptr<ParserType>& reader) {
        while(1) {
            auto state = advance(reader);
            if (state.end) {
                break;
            }
        }

        return _n_reads;
    }

    /**
     * @Synopsis     Default paired-end sequence processing implementation:
     *               just consume both.
     *
     * @Param bundle ReadBundle containing the two sequences.
     */
    void process_sequence(parsing::RecordPair& pair) {
        if (pair.first) {
            derived().process_sequence(pair.first.value());
        }
        if (pair.second) {
            derived().process_sequence(pair.second.value());
        }
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
                          << this->_n_reads
                          << ", exception was "
                          << e.what() << std::endl;
            }
            return {};
        }  catch (parsing::InvalidRead& e) {
            if (_verbose) {
                std::cerr << "WARNING: Bad invalid read encountered at "
                          << this->_n_reads
                          << ", exception was "
                          << e.what() << std::endl;
            }
            return {};
        }
    }

    /**
     * @Synopsis  Consume the next FINE_INTERVAL sequences and return the
     *            current interval state when completed.
     *
     * @Param reader The reader to consume from.
     *
     * @Returns   interval_state with current interval.
     */
    interval_state advance(std::shared_ptr<parsing::SplitPairedReader<ParserType>>& reader) {
        std::optional<parsing::RecordPair> bundle;
        while(!reader->is_complete()) {
            bundle = handle_next(*reader);
            
            if (!bundle) {
                continue;
            }

            derived().process_sequence(bundle.value());
            int _bundle_count = (bool)bundle.value().first + (bool)bundle.value().second;
            _n_reads += _bundle_count;

            auto tick_result = _notify_tick(_bundle_count);
            if (_ticked(tick_result)) {
                return tick_result;
            }
        }
        _notify_stop();
        return interval_state(false, false, false, true);
    }

    /**
     * @Synopsis  Consume the next FINE_INTERVAL sequences and return the
     *            current interval state when completed.
     *
     * @Param parser The reader to consume from.
     *
     * @Returns   interval_state with current interval.
     */
    interval_state advance(std::shared_ptr<ParserType>& parser) {

        std::optional<parsing::Record> record;
        // Iterate through the reads and consume their k-mers.
        int i = 0;
        while (!parser->is_complete()) {
            record = handle_next(*parser);
            
            if(!record) {
                continue;
            }
            
            if(!record) {
                continue;
            }

            derived().process_sequence(record.value());
 
            __sync_add_and_fetch( &_n_reads, 1 );
            auto tick_result = _notify_tick(1);

            if (_ticked(tick_result)) {
                return tick_result;
            }

        }
        _notify_stop();
        return interval_state(false, false, false, true);
    }

    /**
     * @Synopsis  Number of sequences / reads consumed.
     *
     * @Returns   
     */
    uint64_t n_reads() const {
        return _n_reads;
    }

    template<typename... Args>
    static std::shared_ptr<Derived> build(Args&&... args,
                                          uint64_t fine_interval,
                                          uint64_t medium_interval,
                                          uint64_t coarse_interval,
                                          bool     verbose) {
        return Derived::build(std::forward<Args>(args)...,
                              fine_interval,
                              medium_interval,
                              coarse_interval,
                              verbose);
    }

private:

    FileProcessor() : _n_reads(0) {}

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
    uint64_t _n_kmers;

    typedef FileProcessor<InserterProcessor<InserterType, ParserType>,
                          ParserType> Base;

public:

    using Base::process_sequence;
    typedef typename Base::alphabet alphabet;
    
    InserterProcessor(std::shared_ptr<InserterType> inserter,
                      uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                      uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                      uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE,
                      bool     verbose         = false)
        : Base(fine_interval, medium_interval, coarse_interval, verbose),
          inserter(inserter),
          _n_kmers(0)
    {
    }

    void process_sequence(const parsing::Record& read) {
        try {
            inserter->insert_sequence(read.sequence);
        } catch (SequenceLengthException &e) {
            if (this->_verbose) {
                std::cerr << "WARNING: Skipped sequence that was too short: read "
                          << this->_n_reads << " with sequence "
                          << read.sequence 
                          << std::endl;
            }
            return;
        } catch (InvalidCharacterException& e) {
            return;
        } catch (std::exception &e) {
            std::cerr << "ERROR: Exception thrown at " << this->_n_reads 
                      << " with msg: " << e.what()
                      <<  std::endl;
            throw e;
        }
        __sync_add_and_fetch(&_n_kmers, read.sequence.length() - inserter->K + 1);
    }

    void report() {

    }

    uint64_t n_kmers() const {
       return _n_kmers;
    }

    static std::shared_ptr<InserterProcessor<InserterType, ParserType>> build(std::shared_ptr<InserterType> inserter,
                                                                              uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                                                                              uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                                                                              uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE) {
        return std::make_shared<InserterProcessor<InserterType, ParserType>>(inserter, fine_interval, medium_interval, coarse_interval);
    }

};

} //namespace goetia
#endif
