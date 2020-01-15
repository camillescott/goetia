/**
 * (c) Camille Scott, 2019
 * File   : ukhs_signature.cc
 * License: MIT
 * Author : Camille Scott <camille.scott.w@gmail.com>
 * Date   : 31.08.2019
 */


#include "boink/signatures/ukhs_signature.hh"

#include "boink/storage/nibblestorage.hh"
#include "boink/storage/bitstorage.hh"
#include "boink/storage/storage.hh"
#include "boink/storage/qfstorage.hh"
#include "boink/storage/bytestorage.hh"
#include "boink/storage/sparseppstorage.hh"


namespace boink {
namespace signatures {

template <class StorageType, class BaseShifterType>
UnikmerSignature<StorageType, BaseShifterType>::
Reporter::Reporter(std::shared_ptr<UnikmerSignature<StorageType, BaseShifterType>::Signature> signature,
                   const std::string&    filename)
    : SingleFileReporter(filename, "UnikmerSignature::Reporter"),
      signature(signature)
{
    _cerr(this->THREAD_NAME << " reporting at MEDIUM interval.");
    this->msg_type_whitelist.insert(events::MSG_TIME_INTERVAL);
}


template <class StorageType, class BaseShifterType>
void
UnikmerSignature<StorageType, BaseShifterType>::
Reporter::handle_msg(std::shared_ptr<events::Event> event) {
     if (event->msg_type == events::MSG_TIME_INTERVAL) {
        auto _event = static_cast<events::TimeIntervalEvent*>(event.get());
        if (_event->level == events::TimeIntervalEvent::MEDIUM ||
            _event->level == events::TimeIntervalEvent::END) {
            
            _output_stream << _event->t;
            auto counts = signature->get_signature();
            for (auto& count : counts) {
                _output_stream << ", " << count;
            }
            _output_stream << std::endl;
        }
    }       
}

}

template class signatures::UnikmerSignature<storage::BitStorage, hashing::FwdRollingShifter>;
template class signatures::UnikmerSignature<storage::BitStorage, hashing::CanRollingShifter>;
template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::FwdRollingShifter>;
template class signatures::UnikmerSignature<storage::SparseppSetStorage, hashing::CanRollingShifter>;
template class signatures::UnikmerSignature<storage::ByteStorage, hashing::FwdRollingShifter>;
template class signatures::UnikmerSignature<storage::ByteStorage, hashing::CanRollingShifter>;
template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::FwdRollingShifter>;
template class signatures::UnikmerSignature<storage::NibbleStorage, hashing::CanRollingShifter>;
template class signatures::UnikmerSignature<storage::QFStorage, hashing::FwdRollingShifter>;
template class signatures::UnikmerSignature<storage::QFStorage, hashing::CanRollingShifter>;

}

