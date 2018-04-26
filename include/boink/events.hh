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

#ifndef EVENTS_HH
#define EVENTS_HH

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

#include "boink.hh"
#include "event_types.hh"


# ifdef DEBUG_EVENTS
#   define pdebug(x) do { std::ostringstream stream; \
                          stream << std::endl << "@ " << __FILE__ <<\
                          ":" << __FUNCTION__ << ":" <<\
                          __LINE__  << std::endl << x << std::endl;\
                          std::cerr << stream.str(); \
                          } while (0)
# else
#   define pdebug(x) do {} while (0)
# endif

using namespace boink;
using namespace boink::event_types;

namespace boink {
namespace events {

using std::shared_ptr;
using std::make_shared;
using std::unique_ptr;
using std::make_unique;


class ScopedThread {

    std::thread t;

public:

    explicit ScopedThread(std::thread t_)
        : t(std::move(t_)) {

        if(!t.joinable()) {
            throw std::logic_error("Invalid thread.");
        }
    }

    ~ScopedThread() {
        pdebug("Shutting down thread.");
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

    std::unique_ptr<ScopedThread> m_thread;
    std::deque<shared_ptr<Event>> m_queue;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    std::set<event_t> msg_type_whitelist;

public:

    const std::string THREAD_NAME;

    EventListener(const std::string& thread_name)
        : THREAD_NAME(thread_name)
    {
        msg_type_whitelist.insert(MSG_EXIT_THREAD);
        m_thread = make_unique<ScopedThread>
            (std::move(std::thread(&EventListener::process, this)));
    }

    virtual ~EventListener() {
        // puts the poison pill on the queue
        exit_thread();
        // thread joined through RAII
    }

    void exit_thread() {
		auto _exit_event = make_shared<Event>(MSG_EXIT_THREAD);
		notify(_exit_event);
	}

    /// Get the ID of this thread instance
    /// @return The worker thread ID
    std::thread::id get_thread_id() {
        //ASSERT_TRUE(m_thread != 0);
        return m_thread->get_id();
	}

    /// Get the ID of the currently executing thread
    /// @return The current thread ID
    static std::thread::id get_current_thread_id() {
        return std::this_thread::get_id();
	}

    void notify(shared_ptr<Event> event) {
        //std::this_thread::sleep_for(std::chrono::milliseconds(250));
        if (msg_type_whitelist.count(event->msg_type)) {
            std::unique_lock<std::mutex> lk(m_mutex);
            m_queue.push_back(event);
            m_cv.notify_one();
        } else {
            pdebug("Filtered event of type " << event->msg_type);
        }
	}

protected:

    EventListener(const EventListener&);
    EventListener& operator=(const EventListener&);

	void process() {
        pdebug( "EventListener::" << THREAD_NAME << " listening "
                  << "at thread ID " << get_current_thread_id());

		while (1) {
			shared_ptr<Event> msg;
			{
				// Wait for a message to be added to the queue
				std::unique_lock<std::mutex> lk(m_mutex);
				while (m_queue.empty())
					m_cv.wait(lk);

				if (m_queue.empty())
					continue;

				msg = m_queue.front();
				m_queue.pop_front();
			}

			if (msg->msg_type == MSG_EXIT_THREAD) {
				std::unique_lock<std::mutex> lk(m_mutex);
                m_queue.clear();
                pdebug("Exit EventListener::" << THREAD_NAME 
                       << " listener at thread ID " << get_current_thread_id());
                return;

	        }

            handle_msg(msg);
		}
	}

    virtual void handle_msg(shared_ptr<Event> msg) {
        pdebug("Handled no-op message");
    }
};


class EventNotifier {

protected:

    std::set<EventListener*> registered_listeners;

public:

    EventNotifier()
    {
    }

    void notify(shared_ptr<Event> event) {
        //pdebug("Notifying " << registered_listeners.size()
        //       << " listeners of event of type " << event->msg_type);
        for (auto listener : registered_listeners) {
            if (listener != nullptr) {
                listener->notify(event);
            }
        }
    }

    void register_listener(EventListener* listener) {
        pdebug("Register " << listener->THREAD_NAME);
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

#undef pdebug
#endif
