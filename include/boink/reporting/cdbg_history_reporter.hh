/* cdbg_history_reporter.hh -- async reporters
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_CDBG_HISTORY_REPORTER_HH
#define BOINK_CDBG_HISTORY_REPORTER_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "boink/boink.hh"
#include "boink/event_types.hh"
#include "boink/cdbg/cdbg_types.hh"
#include "boink/reporting/reporters.hh"
#include "boink/reporting/report_types.hh"

#include "sparsepp/spp.h"


namespace boink {
namespace reporting {

using boink::cdbg::id_t;

class cDBGHistoryReporter : public SingleFileReporter {
private:

    id_t _edge_id_counter;
    spp::sparse_hash_map<id_t, std::vector<std::string>> node_history;

public:
    cDBGHistoryReporter(const std::string& filename)
        : SingleFileReporter(filename, "cDBGHistoryReporter"),
          _edge_id_counter(0)
    {
        _cerr(this->THREAD_NAME << " reporting continuously.");

        this->msg_type_whitelist.insert(events::MSG_HISTORY_NEW);
        this->msg_type_whitelist.insert(events::MSG_HISTORY_SPLIT);
        this->msg_type_whitelist.insert(events::MSG_HISTORY_SPLIT_CIRCULAR);
        this->msg_type_whitelist.insert(events::MSG_HISTORY_MERGE);
        this->msg_type_whitelist.insert(events::MSG_HISTORY_EXTEND);
        this->msg_type_whitelist.insert(events::MSG_HISTORY_CLIP);
        this->msg_type_whitelist.insert(events::MSG_HISTORY_DELETE);

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

    virtual void handle_exit() {
        _output_stream << "</graph>" << std::endl;
        _output_stream << "</graphml>" << std::endl;
    }

    void write_node(std::string id, id_t boink_id, std::string node_meta, std::string sequence) {
        _output_stream << "<node id=\"" << id << "\">" << std::endl
                       << "    <data key=\"seq\">" << sequence << "</data>" << std::endl
                       << "    <data key=\"meta\">" << node_meta << "</data>" << std::endl
                       << "    <data key=\"node_id\">" << boink_id << "</data>" << std::endl
                       << "</node>" << std::endl;
    }

    void write_edge(std::string src, std::string dst, std::string op) {
        auto id = _edge_id_counter++;
        _output_stream << "<edge id=\"" << id << "\" source=\"" 
                       << src << "\" target=\"" << dst << "\">" << std::endl
                       << "    <data key=\"op\">" << op << "</data>" << std::endl
                       << "</edge>" << std::endl;
    }

    std::string add_node_edit(id_t node_id, cdbg::node_meta_t meta, std::string sequence) {
        auto change_num = node_history[node_id].size();
        std::string id = std::to_string(node_id) + "_" + std::to_string(change_num);
        node_history[node_id].push_back(id);
        write_node(id, node_id, std::string(node_meta_repr(meta)), sequence);
        return id;
    }

    std::string add_new_node(id_t node_id, cdbg::node_meta_t meta, std::string sequence) {
        std::string id = std::to_string(node_id) + "_0";
        if (node_history.count(node_id) == 0) {
            node_history[node_id] = std::vector<std::string>{id};
            write_node(id, node_id, std::string(node_meta_repr(meta)), sequence);
        }
        return id;
    }

    virtual void handle_msg(shared_ptr<events::Event> event) {
        if (event->msg_type == events::MSG_HISTORY_NEW) {
            auto _event = static_cast<events::HistoryNewEvent*>(event.get());
            add_new_node(_event->id, _event->meta, _event->sequence);

        } else if (event->msg_type == events::MSG_HISTORY_SPLIT) {
            auto _event = static_cast<events::HistorySplitEvent*>(event.get());
            
            std::string parent_id = node_history[_event->parent].back();
            std::string lid, rid;
            if (_event->lchild == _event->parent) {
                lid = add_node_edit(_event->lchild, _event->lmeta, _event->lsequence);
                rid = add_new_node(_event->rchild, _event->rmeta, _event->rsequence);
            } else {
                lid = add_new_node(_event->lchild, _event->lmeta, _event->lsequence);
                rid = add_node_edit(_event->rchild, _event->rmeta, _event->rsequence);
            }
            write_edge(parent_id, lid, std::string("SPLIT"));
            write_edge(parent_id, rid, std::string("SPLIT"));

        } else if (event->msg_type == events::MSG_HISTORY_MERGE) {
            auto _event = static_cast<events::HistoryMergeEvent*>(event.get());

            std::string l_parent_id = node_history[_event->lparent].back();
            std::string r_parent_id = node_history[_event->rparent].back();
            std::string child_id = add_node_edit(_event->child, _event->meta, _event->sequence);
            
            write_edge(l_parent_id, child_id, std::string("MERGE"));
            write_edge(r_parent_id, child_id, std::string("MERGE"));

        } else if (event->msg_type == events::MSG_HISTORY_EXTEND) {
            auto _event = static_cast<events::HistoryExtendEvent*>(event.get());

            std::string src = node_history[_event->id].back();
            std::string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
            write_edge(src, dst, std::string("EXTEND"));

        } else if (event->msg_type == events::MSG_HISTORY_CLIP) {
            auto _event = static_cast<events::HistoryClipEvent*>(event.get());

            std::string src = node_history[_event->id].back();
            std::string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
            write_edge(src, dst, std::string("CLIP"));

        } else if (event->msg_type == events::MSG_HISTORY_SPLIT_CIRCULAR) {
            auto _event = static_cast<events::HistorySplitCircularEvent*>(event.get());

            std::string src = node_history[_event->id].back();
            std::string dst = add_node_edit(_event->id, _event->meta, _event->sequence);
            write_edge(src, dst, std::string("SPLIT_CIRCULAR"));
        }
    }
};




}
}

#endif


