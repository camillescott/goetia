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


namespace boink {
namespace events {


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

    std::mutex m_mutex;
    std::condition_variable m_cv;
    std::unique_ptr<ScopedThread> m_thread;
    std::deque<std::shared_ptr<events::Event>> m_queue;
    std::set<events::event_t> msg_type_whitelist;
    bool _shutdown;
    uint64_t _to_process;
    uint64_t MAX_EVENTS;
    uint64_t MIN_EVENTS_RESTART;

public:

    const std::string THREAD_NAME;

    EventListener(const std::string& thread_name)
        : m_mutex(),
          m_cv(),
          _shutdown(false),
          _to_process(0),
          MAX_EVENTS(50000),
          MIN_EVENTS_RESTART(MAX_EVENTS * 0.9),
          THREAD_NAME("EventListener::" + thread_name)
    {
        msg_type_whitelist.insert(MSG_EXIT_THREAD);
        m_thread = std::make_unique<ScopedThread>
            (std::move(std::thread(&EventListener::process, this)));
    }

    virtual ~EventListener() {
        // puts the poison pill on the queue
        if (!_shutdown) {
            exit_thread();
        }
        // thread joined through RAII
    }

    void exit_thread() {
		auto _exit_event = std::make_shared<Event>(MSG_EXIT_THREAD);
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_queue.push_back(_exit_event);
        }
        m_cv.notify_all();

        {
            std::unique_lock<std::mutex> lk(m_mutex);
            m_cv.wait(lk, [this]{ return _shutdown; });
        }
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

    void notify(std::shared_ptr<events::Event> event) {
        //std::this_thread::sleep_for(std::chrono::milliseconds(250));
        if (msg_type_whitelist.count(event->msg_type)) {
            {
                std::unique_lock<std::mutex> lk(m_mutex);
                if (_to_process > MAX_EVENTS) {
                    _cerr(THREAD_NAME << " hit MAX_EVENTS ("
                          << MAX_EVENTS << ") on queue; blocking "
                          "until " << MIN_EVENTS_RESTART << " events.");
                    wait_on_processing(MIN_EVENTS_RESTART);

                }
                ++_to_process;
                m_queue.push_back(event);
            }
            m_cv.notify_one();
        } else {
            pdebug("Filtered event of type " << event->msg_type);
        }
	}

    void clear_events() {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_queue.clear();
    }

    void wait_on_processing(uint64_t min_events_restart) {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cv.wait(lock, [&]{ return _to_process <= min_events_restart; });
    }

protected:

    EventListener(const EventListener&);
    EventListener& operator=(const EventListener&);

	void process() {
        _cerr(THREAD_NAME << " listening "
              << "at thread ID " << get_current_thread_id());

		while (1) {
            std::shared_ptr<Event> msg;
			{
				// Wait for a message to be added to the queue
				std::unique_lock<std::mutex> lk(m_mutex);
				m_cv.wait(lk, [this]{ return !m_queue.empty(); });

				if (m_queue.empty())
					continue;

				msg = m_queue.front();
				m_queue.pop_front();

                if (m_queue.size() < MIN_EVENTS_RESTART) {
                    m_cv.notify_one();
                }
			}

			if (msg->msg_type == MSG_EXIT_THREAD) {
                handle_exit();
                _cerr("Exit " << THREAD_NAME 
                      << " listener at thread ID " << get_current_thread_id());
                {
                    std::unique_lock<std::mutex> lk(m_mutex);
                    _shutdown = true;
                }
                m_cv.notify_all();
                return;
	        }

            handle_msg(msg);

            {
                std::unique_lock<std::mutex> lk(m_mutex);
                --_to_process;
                m_cv.notify_all();
            }
		}
	}

    virtual void handle_msg(std::shared_ptr<events::Event> msg) {
        pdebug("Handled no-op message");
    }

    virtual void handle_exit() {
        pdebug("Handled no-op exit.");
    }
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
        _cerr("Register " << listener->THREAD_NAME << " at " << this);
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
