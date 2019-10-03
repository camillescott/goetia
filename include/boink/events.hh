/* events.hh -- event listeners and senders
 * 
 * Code in this file modified and extended from:
 * `https://www.codeproject.com/Articles/1169105/Cplusplus-
 * std-thread-Event-Loop-with-Message-Queue`
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#ifndef BOINK_EVENTS_HH
#define BOINK_EVENTS_HH

#include <iostream>
#include <sstream>
#include <deque>
#include <thread>
#include <memory>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <set>
#include <chrono>

#include "boink/event_types.hh"


namespace boink {
namespace events {

// RAII std::thread wrapper
class ScopedThread {

    std::thread t;

public:

    explicit ScopedThread(std::thread t_);

    ~ScopedThread() {
        t.join();
    }

    std::thread::id get_id() {
        return t.get_id();
    }

    ScopedThread(ScopedThread const&)=delete;
    ScopedThread& operator=(ScopedThread const&)=delete;
};


class EventListener  {

protected:

    std::mutex                                 msg_mutex;
    std::condition_variable                    msg_cv;
    std::unique_ptr<ScopedThread>              listener_thread;
    std::deque<std::shared_ptr<events::Event>> msg_queue;
    std::set<events::event_t>                  msg_type_whitelist;
    bool                                       _shutdown;
    uint64_t                                   _to_process;
    uint64_t                                   MAX_EVENTS;
    uint64_t                                   MIN_EVENTS_RESTART;

public:

    const std::string THREAD_NAME;

    EventListener(const std::string& thread_name);

    virtual ~EventListener();

    void exit_thread();

    /// Get the ID of this thread instance
    /// @return The worker thread ID
    std::thread::id get_thread_id();

    /// Get the ID of the currently executing thread
    /// @return The current thread ID
    static std::thread::id get_current_thread_id();

    void notify(std::shared_ptr<events::Event> event);

    void clear_events();

    void wait_on_processing(uint64_t min_events_restart);

protected:

    //EventListener(const EventListener&);
    //EventListener& operator=(const EventListener&);

	void process();

    virtual void handle_msg(std::shared_ptr<events::Event> msg);

    virtual void handle_exit();
};


class EventNotifier {

protected:

    std::set<EventListener*> registered_listeners;

public:

    EventNotifier()
    {
    }

    ~EventNotifier()
    {
        stop_listeners();
    }

    void notify(std::shared_ptr<events::Event> event) {
        //pdebug("Notifying " << registered_listeners.size()
        //       << " listeners of event of type " << event->msg_type);
        for (auto listener : registered_listeners) {
            if (listener != nullptr) {
                listener->notify(event);
            }
        }
    }

    void register_listener(EventListener* listener) {
        //_cerr("Register " << listener->THREAD_NAME << " at " << this);
        registered_listeners.insert(listener);
    }

    void stop_listeners() {
        for (auto listener : registered_listeners) {
            if (listener != nullptr) {
                listener->exit_thread();
            }
        }
        registered_listeners.clear();
    }

    void clear_listeners() {
        registered_listeners.clear();
    }

};



}
}

#endif
