/* reporters.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef REPORTERS_HH
#define REPORTERS_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "boink/boink.hh"
#include "boink/cdbg.hh"
#include "boink/compactor.hh"
#include "boink/events.hh"
#include "boink/event_types.hh"

#include <sparsepp/sparsepp/spp.h>

namespace boink {
namespace reporters {

using std::shared_ptr;

using boink::events::EventListener;
using boink::event_types::Event;
using boink::event_types::TimeIntervalEvent;

class SingleFileReporter : public EventListener {

protected:

    std::string _output_filename;
    std::ofstream _output_stream;

public:

    SingleFileReporter(const std::string& output_filename,
             const std::string& thread_name) 
        : EventListener(thread_name),
          _output_filename(output_filename),
          _output_stream(_output_filename.c_str())
    {
    }

    SingleFileReporter(const std::string& output_filename)
        : SingleFileReporter(output_filename, "SingleFileReporter")
    {
    }

    virtual ~SingleFileReporter() {
        _output_stream.close();
    }
};


class MultiFileReporter : public EventListener {

protected:

    std::string _file_prefix;
    std::vector<std::string> _filenames;
    std::vector<std::ofstream> _streams;

public:

    MultiFileReporter(const std::string& prefix,
                      const std::string& thread_name)
        : EventListener(thread_name),
          _file_prefix(prefix)
    {
    }

    MultiFileReporter(const std::string& prefix)
        : MultiFileReporter(prefix, "MultiFileReporter")
    {
    }

    std::ofstream& current_stream() {
        return _streams.back();
    }

    std::string& current_filename() {
        return _filenames.back();
    }

    std::ofstream& next_stream(uint64_t start_time, 
                               const std::string& suffix) {
        std::ostringstream name;
        name << _file_prefix << "." << start_time << "."
             << suffix;
        std::string filename = name.str();

        if (_streams.size() > 0) {
            current_stream().close();   
        }
        _streams.push_back(std::move(std::ofstream(filename.c_str())));
        _filenames.push_back(filename);
        return _streams.back();
    }

    virtual ~MultiFileReporter() {
        current_stream().close();
    }

};


template <class GraphType>
class StreamingCompactorReporter: public SingleFileReporter {

protected:

    StreamingCompactor<GraphType> * compactor;

public:

    StreamingCompactorReporter(StreamingCompactor<GraphType> * compactor,
                               const std::string& output_filename)
        : SingleFileReporter(output_filename, "StreamingCompactorReporter"),
          compactor(compactor)
    {    
        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);

        _output_stream << "read_n,n_full,n_tips,n_islands,n_trivial"
                          ",n_circular,n_loops,n_dnodes,n_unodes,n_tags,"
                          "n_updates,n_unique,estimated_fp" << std::endl;
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::FINE ||
                _event->level == TimeIntervalEvent::END) {
                auto report = compactor->get_report();
                _output_stream << _event->t << ","
                               << report->n_full << ","
                               << report->n_tips << ","
                               << report->n_islands << ","
                               << report->n_trivial << ","
                               << report->n_circular << ","
                               << report->n_loops << ","
                               << report->n_dnodes << ","
                               << report->n_unodes << ","
                               << report->n_tags << ","
                               << report->n_updates << ","
                               << report->n_unique << ","
                               << report->estimated_fp 
                               << std::endl;
            }
        }
    }
};


class cDBGHistoryReporter : public SingleFileReporter {
private:

    id_t _edge_id_counter;
    spp::sparse_hash_map<id_t, std::vector<string>> node_history;

public:
    cDBGHistoryReporter(const std::string& filename)
        : SingleFileReporter(filename, "cDBGHistoryReporter"),
          _edge_id_counter(0)
    {
        this->msg_type_whitelist.insert(boink::event_types::MSG_DAG_NEW);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DAG_SPLIT);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DAG_SPLIT_CIRCULAR);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DAG_MERGE);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DAG_EXTEND);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DAG_CLIP);
        this->msg_type_whitelist.insert(boink::event_types::MSG_DAG_DELETE);

        _output_stream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                          "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" "
                          "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                          "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns "
                          "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">"
                       << std::endl // the header, open <graphml>
                       << "<graph id=\"cDBG_History_DAG\" edgedefault=\"directed\">" << std::endl
                       << "<key id=\"op\" for=\"edge\" attr.name=\"op\" attr.type=\"string\"/>" << std::endl
                       << "<key id=\"seq\" for=\"node\" attr.name=\"seq\" attr.type=\"string\"/>" << std::endl
                       << "<key id=\"meta\" for=\"node\" attr.name=\"meta\" attr.type=\"string\"/>" << std::endl
                       << "<key id=\"node_id\" for=\"node\" attr.name=\"node_id\" attr.type=\"long\"/>"
                       << std::endl; // open <graph>
    }

    virtual ~cDBGHistoryReporter() {
        _output_stream << "</graph>" << std::endl;
        _output_stream << "</graphml>" << std::endl;
    }

    void write_node(string id, id_t boink_id, string node_meta, string sequence) {
        _output_stream << "<node id=\"" << id << "\">" << std::endl
                       << "    <data key=\"seq\">" << sequence << "</data>" << std::endl
                       << "    <data key=\"meta\">" << node_meta << "</data>" << std::endl
                       << "    <data key=\"node_id\">" << boink_id << "</data>" << std::endl
                       << "</node>" << std::endl;
    }

    void write_edge(string src, string dst, string op) {
        auto id = _edge_id_counter++;
        _output_stream << "<edge id=\"" << id << "\" source=\"" 
                       << src << "\" target=\"" << dst << "\">" << std::endl
                       << "    <data key=\"op\">" << op << "</data>" << std::endl
                       << "</edge>" << std::endl;
    }

    string add_node_edit(id_t node_id, node_meta_t meta, string sequence) {
        auto change_num = node_history[node_id].size();
        string id = std::to_string(node_id) + "_" + std::to_string(change_num);
        node_history[node_id].push_back(id);
        write_node(id, node_id, string(node_meta_repr(meta)), sequence);
        return id;
    }

    string add_new_node(id_t node_id, node_meta_t meta, string sequence) {
        string id = std::to_string(node_id) + "_0";
        if (node_history.count(node_id) == 0) {
            node_history[node_id] = std::vector<string>{id};
            write_node(id, node_id, string(node_meta_repr(meta)), sequence);
        }
        return id;
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_DAG_NEW) {
            auto _event = static_cast<DAGNewEvent*>(event.get());
            add_new_node(_event->id, _event->meta, _event->sequence);

        } else if (event->msg_type == boink::event_types::MSG_DAG_SPLIT) {
            auto _event = static_cast<DAGSplitEvent*>(event.get());
            
            string parent_id = node_history[_event->parent].back();
            string lid, rid;
            if (_event->lchild == _event->parent) {
                lid = add_node_edit(_event->lchild, _event->lmeta, _event->lsequence);
                rid = add_new_node(_event->rchild, _event->rmeta, _event->rsequence);
            } else {
                lid = add_new_node(_event->lchild, _event->lmeta, _event->lsequence);
                rid = add_node_edit(_event->rchild, _event->rmeta, _event->rsequence);
            }
            write_edge(parent_id, lid, string("SPLIT"));
            write_edge(parent_id, rid, string("SPLIT"));

        } else if (event->msg_type == boink::event_types::MSG_DAG_MERGE) {
            auto _event = static_cast<DAGMergeEvent*>(event.get());

            string l_parent_id = node_history[_event->lparent].back();
            string r_parent_id = node_history[_event->rparent].back();
            string child_id = add_node_edit(_event->child, _event->meta, _event->sequence);
            
            write_edge(l_parent_id, child_id, string("MERGE"));
            write_edge(r_parent_id, child_id, string("MERGE"));

        } else if (event->msg_type == boink::event_types::MSG_DAG_EXTEND) {
            auto _event = static_cast<DAGExtendEvent*>(event.get());

            string src = node_history[_event->id].back();
            string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
            write_edge(src, dst, string("EXTEND"));

        } else if (event->msg_type == boink::event_types::MSG_DAG_CLIP) {
            auto _event = static_cast<DAGClipEvent*>(event.get());

            string src = node_history[_event->id].back();
            string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
            write_edge(src, dst, string("CLIP"));

        } else if (event->msg_type == boink::event_types::MSG_DAG_SPLIT_CIRCULAR) {
            auto _event = static_cast<DAGSplitCircularEvent*>(event.get());

            string src = node_history[_event->id].back();
            string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
            write_edge(src, dst, string("SPLIT_CIRCULAR"));
        }
    }
};


template <class GraphType>
class cDBGWriter : public MultiFileReporter {
protected:

    cDBG<GraphType> * cdbg;
    cDBGFormat format;

public:

    cDBGWriter(cDBG<GraphType> * cdbg,
               cDBGFormat format,
               const string& output_prefix)
        : MultiFileReporter(output_prefix,
                            "cDBGWriter[" + cdbg_format_repr(format) + "]"),
          cdbg(cdbg),
          format(format)
    {
        this->msg_type_whitelist.insert(boink::event_types::MSG_TIME_INTERVAL);
    }

    virtual void handle_msg(shared_ptr<Event> event) {
        if (event->msg_type == boink::event_types::MSG_TIME_INTERVAL) {
            auto _event = static_cast<TimeIntervalEvent*>(event.get());
            if (_event->level == TimeIntervalEvent::MEDIUM ||
                _event->level == TimeIntervalEvent::END) {

                std::ofstream& stream = this->next_stream(_event->t,
                                                          cdbg_format_repr(format));
                std::string&   filename = this->current_filename();

                _cerr(this->THREAD_NAME << ", t=" << _event->t <<
                      ": write cDBG to " << filename);
                cdbg->write(stream, format);
            }
        }
    }

};

}
}

#endif
