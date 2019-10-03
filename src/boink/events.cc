/* events.cc
 *
 * Copyright (C) 2018 Camille Scott
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

#include <thread>

#include "boink/events.hh"

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


ScopedThread::ScopedThread(std::thread t_) 
    : t(std::move(t_))
{
    
    if(!t.joinable()) {
        throw std::logic_error("Invalid thread.");
    }
}

EventListener::~EventListener() {
    // puts the poison pill on the queue
    if (!_shutdown) {
        exit_thread();
    }
    // thread joined through RAII
}


EventListener::EventListener(const std::string& thread_name)
    : msg_mutex(),
      msg_cv(),
      _shutdown(false),
      _to_process(0),
      MAX_EVENTS(50000),
      MIN_EVENTS_RESTART(MAX_EVENTS * 0.9),
      THREAD_NAME("EventListener::" + thread_name)
{
    msg_type_whitelist.insert(MSG_EXIT_THREAD);
    listener_thread = std::make_unique<ScopedThread>
        (std::move(std::thread(&EventListener::process, this)));
}


void EventListener::exit_thread() {
    auto _exit_event = std::make_shared<Event>(MSG_EXIT_THREAD);
    {
        std::unique_lock<std::mutex> lock(msg_mutex);
        msg_queue.push_back(_exit_event);
    }
    msg_cv.notify_all();

    {
        std::unique_lock<std::mutex> lk(msg_mutex);
        msg_cv.wait(lk, [this]{ return _shutdown; });
    }
}


std::thread::id EventListener::get_thread_id() {
    //ASSERT_TRUE(listener_thread != 0);
    return listener_thread->get_id();
}

std::thread::id EventListener::get_current_thread_id() {
    return std::this_thread::get_id();
}

void EventListener::notify(std::shared_ptr<events::Event> event) {
    //std::this_thread::sleep_for(std::chrono::milliseconds(250));
    if (msg_type_whitelist.count(event->msg_type)) {
        {
            std::unique_lock<std::mutex> lk(msg_mutex);
            if (_to_process > MAX_EVENTS) {
                _cerr(THREAD_NAME << " hit MAX_EVENTS ("
                      << MAX_EVENTS << ") on queue; blocking "
                      "until " << MIN_EVENTS_RESTART << " events.");
                wait_on_processing(MIN_EVENTS_RESTART);

            }
            ++_to_process;
            msg_queue.push_back(event);
        }
        msg_cv.notify_one();
    } else {
        pdebug("Filtered event of type " << event->msg_type);
    }
}


void EventListener::clear_events() {
    std::unique_lock<std::mutex> lock(msg_mutex);
    msg_queue.clear();
}

void EventListener::wait_on_processing(uint64_t min_events_restart) {
    std::unique_lock<std::mutex> lock(msg_mutex);
    msg_cv.wait(lock, [&]{ return _to_process <= min_events_restart; });
}

void EventListener::process() {
    _cerr(THREAD_NAME << " listening "
          << "at thread ID " << get_current_thread_id());

    while (1) {
        std::shared_ptr<Event> msg;
        {
            // Wait for a message to be added to the queue
            std::unique_lock<std::mutex> lk(msg_mutex);
            msg_cv.wait(lk, [this]{ return !msg_queue.empty(); });

            if (msg_queue.empty())
                continue;

            msg = msg_queue.front();
            msg_queue.pop_front();

            if (msg_queue.size() < MIN_EVENTS_RESTART) {
                msg_cv.notify_one();
            }
        }

        if (msg->msg_type == MSG_EXIT_THREAD) {
            handle_exit();
            _cerr("Exit " << THREAD_NAME 
                  << " listener at thread ID " << get_current_thread_id());
            {
                std::unique_lock<std::mutex> lk(msg_mutex);
                _shutdown = true;
            }
            msg_cv.notify_all();
            return;
        }

        handle_msg(msg);

        {
            std::unique_lock<std::mutex> lk(msg_mutex);
            --_to_process;
            msg_cv.notify_all();
        }
    }
}

void EventListener::handle_msg(std::shared_ptr<events::Event> msg) {
    pdebug("Handled no-op message");
}

void EventListener::handle_exit() {
    pdebug("Handled no-op exit.");
}


}
}

#undef pdebug
