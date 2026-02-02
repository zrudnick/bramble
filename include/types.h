
#pragma once

#include <ankerl/unordered_dense.h>

using tid_t = uint32_t;
using pos_t = uint32_t;
using read_id_t = uint32_t;
using refid_t = int32_t;

// Use ankerl::unordered_dense for better performance
template<typename Key, typename Value>
using unordered_map = ankerl::unordered_dense::map<Key, Value>;

template<typename Key>
using unordered_set = ankerl::unordered_dense::set<Key>;