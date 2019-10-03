/* processors.hh -- 
#include "boink/hashing/hashing_types.hh"
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_PROCESSORS_HH
#define BOINK_PROCESSORS_HH

#include <tuple>
#include <memory>
#include <fstream>
#include <iostream>
#include <string>

#include "boink/boink.hh"
#include "boink/parsing/parsing.hh"
#include "boink/parsing/readers.hh"
#include "boink/events.hh"
#include "boink/event_types.hh"
#include "boink/hashing/exceptions.hh"


class KmerMinHash;

namespace boink {

struct DEFAULT_INTERVALS {

    constexpr static unsigned long FINE = 10000;
    constexpr static unsigned long MEDIUM = 100000;
    constexpr static unsigned long COARSE = 1000000;

};


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
        if (counter == interval) {
            counter = 0;
            return true;
        } else {
            return false;
        }
    }

};


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


template <class Derived,
          class ParserType = parsing::FastxReader>
class FileProcessor : public events::EventNotifier {

protected:

    std::array<IntervalCounter, 3> counters;
    uint64_t _n_reads;

    bool _ticked(interval_state tick) {
        return tick.fine || tick.medium || tick.coarse || tick.end;
    }

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

    void _notify_stop() {
        auto event = std::make_shared<events::TimeIntervalEvent>();
        event->level = events::TimeIntervalEvent::END;
        event->t = _n_reads;
        notify(event);
    }


public:

    using events::EventNotifier::register_listener;
    using events::EventNotifier::notify;

    FileProcessor(uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                  uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                  uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE)
        :  events::EventNotifier(),
           counters ({{ fine_interval, 
                        medium_interval,
                        coarse_interval }}),
          _n_reads(0) {
        
        //std::cout << counters[0].counter << " " << counters[1].counter << std::endl;
    }

    uint64_t process(const std::string& left_filename,
                     const std::string& right_filename,
                     uint32_t min_length=0,
                     bool force_name_match=false) {
        parsing::SplitPairedReader<ParserType> reader(left_filename,
                                             right_filename,
                                             min_length,
                                             force_name_match);
        return process(reader);
    }

    uint64_t process(std::string const &filename) {
        parsing::ReadParserPtr<ParserType> parser = parsing::get_parser<ParserType>(filename);
        return process(parser);
    }

    uint64_t process(parsing::SplitPairedReader<ParserType>& reader) {
        while(1) {
            auto state = advance(reader);
            if (state.end) {
                break;
            }
        }

        return _n_reads;
    }

    uint64_t process(parsing::ReadParserPtr<ParserType>& parser) {
        while(1) {
            auto state = advance(parser);
            if (state.end) {
                break;
            }
        }

        return _n_reads;
    }

    void process_sequence(parsing::ReadBundle& bundle) {
        if (bundle.has_left) {
            derived().process_sequence(bundle.left);
        }
        if (bundle.has_right) {
            derived().process_sequence(bundle.right);
        }
    }

    interval_state advance(parsing::SplitPairedReader<ParserType>& reader) {
        parsing::ReadBundle bundle;

        while(!reader.is_complete()) {
            bundle = reader.next();
            derived().process_sequence(bundle);

            int _bundle_count = bundle.has_left + bundle.has_right;
            _n_reads += _bundle_count;

            auto tick_result = _notify_tick(_bundle_count);
            if (_ticked(tick_result)) {
                return tick_result;
            }
        }
        _notify_stop();
        return interval_state(false, false, false, true);
    }

    interval_state advance(parsing::ReadParserPtr<ParserType>& parser) {
        parsing::Read read;

        // Iterate through the reads and consume their k-mers.
        while (!parser->is_complete()) {
            try {
                read = parser->get_next_read( );
            } catch (parsing::NoMoreReadsAvailable) {
                break;
            }

            read.set_clean_seq();
            derived().process_sequence(read);
 

            __sync_add_and_fetch( &_n_reads, 1 );
            auto tick_result = _notify_tick(1);

            if (_ticked(tick_result)) {
                return tick_result;
            }

        }
        _notify_stop();
        return interval_state(false, false, false, true);
    }

    uint64_t n_reads() const {
        return _n_reads;
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


template <class InserterType,
          class ParserType = parsing::FastxReader>
class InserterProcessor : public FileProcessor<InserterProcessor<InserterType, ParserType>,
                                               ParserType> {

protected:

    std::shared_ptr<InserterType> inserter;
    uint64_t _n_inserted;

    typedef FileProcessor<InserterProcessor<InserterType, ParserType>,
                          ParserType> Base;

public:

    using Base::process_sequence;
    
    InserterProcessor(std::shared_ptr<InserterType> inserter,
                      uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                      uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                      uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE)
        : Base(fine_interval, medium_interval, coarse_interval),
          inserter(inserter), _n_inserted(0) {

    }

    void process_sequence(const parsing::Read& read) {
        auto this_n_inserted = inserter->insert_sequence(read.cleaned_seq);
        __sync_add_and_fetch(&_n_inserted, this_n_inserted);
    }

    void report() {
    }

    uint64_t n_inserted() const {
        return _n_inserted;
    }

    static std::shared_ptr<InserterProcessor<InserterType, ParserType>> build(std::shared_ptr<InserterType> inserter,
                                                                              uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                                                                              uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                                                                              uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE) {
        return std::make_shared<InserterProcessor<InserterType, ParserType>>(inserter, fine_interval, medium_interval, coarse_interval);
    }

};



class SourmashSignatureProcessor : public FileProcessor<SourmashSignatureProcessor,
                                                        parsing::FastxReader> {

protected:

    KmerMinHash * signature;
    typedef FileProcessor<SourmashSignatureProcessor, parsing::FastxReader> Base;

public:

    using Base::process_sequence;

    SourmashSignatureProcessor(KmerMinHash * signature,
                               uint64_t fine_interval   = DEFAULT_INTERVALS::FINE,
                               uint64_t medium_interval = DEFAULT_INTERVALS::MEDIUM,
                               uint64_t coarse_interval = DEFAULT_INTERVALS::COARSE);

    void process_sequence(const parsing::Read& read);

    void report();
};


} //namespace boink
#endif
