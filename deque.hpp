#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <vector>

const unsigned int kBlockSize = 1000;
template <typename T, typename Allocator = std::allocator<T>>
class Deque {
 private:
  std::vector<T*> elems_{};
  size_t first_pos_ = 0;
  size_t last_pos_ = 0;
  Allocator alloc_ = Allocator();
  using alloc = typename std::allocator_traits<Allocator>;

 public:
  template <bool IsConst>
  class common_iterator;

  Deque() {
    std::vector<T*> tmp;
    int idx = 0;
    try {
      tmp.resize(3);
      for (; idx < 3; ++idx) {
        tmp[idx] = alloc::allocate(alloc_, kBlockSize);
      }
    } catch (...) {
      for (int j = 0; j < 3; j++) {
        alloc::deallocate(alloc_, tmp[j], kBlockSize);
      }
      throw;
    }
    first_pos_ = kBlockSize * tmp.size() / 2 - 1;
    last_pos_ = first_pos_;
    std::swap(elems_, tmp);
  }

  Deque(const Deque& to_copy)
      : first_pos_(to_copy.first_pos_), last_pos_(to_copy.last_pos_) {
    alloc_ = alloc::select_on_container_copy_construction(to_copy.alloc_);
    size_t size = elems_.size();
    elems_.resize(to_copy.elems_.size());
    size_t cnt = 0;
    for (size_t i = 0; i < elems_.size(); i++) {
      elems_[i] = alloc::allocate(alloc_, kBlockSize);
      try {
        for (size_t j = 0; j < kBlockSize; j++) {
          if (i * kBlockSize + j >= first_pos_ &&
              i * kBlockSize + j < last_pos_) {
            alloc::construct(alloc_, &elems_[i][j], to_copy.elems_[i][j]);
            cnt++;
          }
        }
      } catch (...) {
        for (size_t j = 0; j < elems_.size(); j++) {
          for (size_t k = 0; k < kBlockSize; k++) {
            if (j * kBlockSize + k >= first_pos_ &&
                j * kBlockSize + k < last_pos_ && cnt != 0) {
              alloc::destroy(alloc_, &elems_[j][k]);
              cnt--;
            }
          }
          alloc::deallocate(alloc_, elems_[j], kBlockSize);
        }
        elems_.resize(size);
        throw;
      }
    }
  }

  Deque(size_t count, const Allocator& alloc = Allocator()) : alloc_(alloc) {
    size_t tmp_size =
        count / kBlockSize + static_cast<size_t>(count % kBlockSize != 0);
    elems_.resize(tmp_size);
    first_pos_ = (kBlockSize * elems_.size() / 2 - 1) - count / 2 -
                 static_cast<size_t>(count % 2 != 0);
    last_pos_ = (kBlockSize * elems_.size() / 2 - 1) + count / 2;
    size_t cnt = 0;
    for (size_t i = 0; i < elems_.size(); ++i) {
      elems_[i] = alloc::allocate(alloc_, kBlockSize);
      try {
        for (size_t j = 0; j < kBlockSize; j++) {
          if (i * kBlockSize + j >= first_pos_ &&
              i * kBlockSize + j < last_pos_) {
            alloc::construct(alloc_, &elems_[i][j]);
            cnt++;
          }
        }
      } catch (...) {
        for (size_t j = 0; j < elems_.size(); j++) {
          for (size_t k = 0; k < kBlockSize; k++) {
            if (j * kBlockSize + k >= first_pos_ &&
                j * kBlockSize + k < last_pos_ && cnt != 0) {
              alloc::destroy(alloc_, &elems_[j][k]);
              cnt--;
            }
          }
          alloc::deallocate(alloc_, elems_[j], kBlockSize);
        }
        first_pos_ = last_pos_ = 0;
        throw;
      }
    }
  }

  Deque(size_t count, const T& value, const Allocator& alloc = Allocator())
      : alloc_(alloc) {
    size_t tmp_size =
        count / kBlockSize * 2 + static_cast<size_t>(count % kBlockSize != 0);
    first_pos_ = (kBlockSize * tmp_size / 2 - 1) - count / 2 -
                 static_cast<size_t>(count % 2 != 0);
    last_pos_ = (kBlockSize * tmp_size / 2 - 1) + count / 2;
    elems_.resize(tmp_size);
    size_t cnt = 0;
    for (size_t i = 0; i < elems_.size(); ++i) {
      elems_[i] = alloc::allocate(alloc_, kBlockSize);
      try {
        for (size_t j = 0; j < kBlockSize; j++) {
          if (i * kBlockSize + j >= first_pos_ &&
              i * kBlockSize + j < last_pos_) {
            alloc::construct(alloc_, &elems_[i][j], value);
            cnt++;
          }
        }
      } catch (...) {
        for (size_t j = 0; j < elems_.size(); j++) {
          for (size_t k = 0; k < kBlockSize; k++) {
            if (j * kBlockSize + k >= first_pos_ &&
                j * kBlockSize + k < last_pos_ && cnt != 0) {
              alloc::destroy(alloc_, &elems_[j][k]);
              cnt--;
            }
          }
          alloc::deallocate(alloc_, elems_[j], kBlockSize);
        }
        first_pos_ = last_pos_ = 0;
        elems_.resize(0);
        throw;
      }
    }
  }

  Deque(Deque&& other)
      : elems_{std::move(other.elems_)},
        first_pos_{std::move(other.first_pos_)},
        last_pos_{std::move(other.last_pos_)},
        alloc_{std::move(other.alloc_)} {
    other.first_pos_ = other.last_pos_ = 0;
  }

  Deque(std::initializer_list<T> init, const Allocator& alloc = Allocator())
      : alloc_(alloc) {
    size_t tmp_size = init.size() / kBlockSize * 2 +
                      static_cast<size_t>(init.size() % kBlockSize != 0);
    first_pos_ = (kBlockSize * tmp_size / 2 - 1) - init.size() / 2 -
                 static_cast<size_t>(init.size() % 2 != 0);
    last_pos_ = (kBlockSize * tmp_size / 2 - 1) + init.size() / 2;
    elems_.resize(tmp_size);
    size_t cnt = 0;
    auto iter = init.begin();
    for (size_t i = 0; i < elems_.size(); ++i) {
      elems_[i] = alloc::allocate(alloc_, kBlockSize);
      try {
        for (size_t j = 0; j < kBlockSize; j++) {
          if (i * kBlockSize + j >= first_pos_ &&
              i * kBlockSize + j < last_pos_) {
            alloc::construct(alloc_, &elems_[i][j], *iter);
            cnt++;
            iter++;
          }
        }
      } catch (...) {
        for (size_t j = 0; j < elems_.size(); j++) {
          for (size_t k = 0; k < kBlockSize; k++) {
            if (j * kBlockSize + k >= first_pos_ &&
                j * kBlockSize + k < last_pos_ && cnt != 0) {
              alloc::destroy(alloc_, &elems_[j][k]);
              cnt--;
            }
          }
          alloc::deallocate(alloc_, elems_[j], kBlockSize);
        }
        first_pos_ = last_pos_ = 0;
        elems_.resize(0);
        throw;
      }
    }
  }

  ~Deque() {
    for (size_t j = 0; j < elems_.size(); j++) {
      for (size_t k = 0; k < kBlockSize; k++) {
        if (j * kBlockSize + k >= first_pos_ &&
            j * kBlockSize + k < last_pos_) {
          alloc::destroy(alloc_, &elems_[j][k]);
        }
      }
      alloc::deallocate(alloc_, elems_[j], kBlockSize);
    }
  }

  Deque& operator=(const Deque& other) {
    if (this != &other) {
      if constexpr (std::allocator_traits<Allocator>::
                        propagate_on_container_copy_assignment::value) {
        alloc_ = other.alloc_;
      }
      Deque tmp(other);
      std::swap(first_pos_, tmp.first_pos_);
      std::swap(elems_, tmp.elems_);
      std::swap(last_pos_, tmp.last_pos_);
    }
    return *this;
  }

  Deque& operator=(Deque&& other) {
    if (this != &other) {
      if constexpr (std::allocator_traits<Allocator>::
                        propagate_on_container_move_assignment::value) {
        alloc_ = std::move(other.alloc_);
      }
      elems_ = std::move(other.elems_);
      first_pos_ = std::move(other.first_pos_);
      last_pos_ = std::move(other.last_pos_);
      other.first_pos_ = other.last_pos_ = 0;
    }
    return *this;
  }

  Allocator get_allocator() const { return alloc_; }

  size_t size() const { return last_pos_ - first_pos_; }

  bool empty() const { return (last_pos_ == first_pos_); }

  T& operator[](size_t idx) {
    size_t pos = first_pos_ + idx;
    return elems_[pos / kBlockSize][pos % kBlockSize];
  }

  const T& operator[](size_t idx) const {
    size_t pos = first_pos_ + idx;
    return elems_[pos / kBlockSize][pos % kBlockSize];
  }

  T& at(size_t idx) {
    size_t pos = first_pos_ + idx;
    if (pos >= last_pos_ || pos < first_pos_) {
      throw std::out_of_range("");
    }
    return elems_[pos / kBlockSize][pos % kBlockSize];
  }

  const T& at(size_t idx) const {
    size_t pos = first_pos_ + idx;
    if (pos >= last_pos_ || pos < first_pos_) {
      throw std::out_of_range("OOE");
    }
    return elems_[pos / kBlockSize][pos % kBlockSize];
  }

  void push_back(const T& value) {
    if (last_pos_ < kBlockSize * elems_.size()) {
      try {
        alloc::construct(
            alloc_, &elems_[last_pos_ / kBlockSize][last_pos_ % kBlockSize],
            value);
      } catch (...) {
        throw;
      }
      last_pos_++;
    } else {
      size_t prev_size = elems_.size();
      elems_.resize(prev_size * 2);
      for (size_t i = prev_size; i < elems_.size(); ++i) {
        elems_[i] = alloc::allocate(alloc_, kBlockSize);
      }
      try {
        alloc::construct(
            alloc_, &elems_[last_pos_ / kBlockSize][last_pos_ % kBlockSize],
            value);
      } catch (...) {
        elems_.resize(prev_size);
        throw;
      }
      last_pos_++;
    }
  }

  void push_back(T&& value) {
    if (last_pos_ < kBlockSize * elems_.size()) {
      try {
        alloc::construct(
            alloc_, &elems_[last_pos_ / kBlockSize][last_pos_ % kBlockSize],
            std::move(value));
      } catch (...) {
        throw;
      }
      last_pos_++;
    } else {
      size_t prev_size = elems_.size();
      elems_.resize(prev_size * 2);
      for (size_t i = prev_size; i < elems_.size(); ++i) {
        elems_[i] = alloc::allocate(alloc_, kBlockSize);
      }
      try {
        alloc::construct(
            alloc_, &elems_[last_pos_ / kBlockSize][last_pos_ % kBlockSize],
            std::move(value));
      } catch (...) {
        elems_.resize(prev_size);
        throw;
      }
      last_pos_++;
    }
  }

  void push_front(const T& value) {
    if (first_pos_ >= 1) {
      first_pos_--;
      try {
        alloc::construct(
            alloc_, &elems_[first_pos_ / kBlockSize][first_pos_ % kBlockSize],
            value);
      } catch (...) {
        throw;
      }
    } else {
      size_t prev_size = elems_.size();
      std::vector<T*> new_elems(prev_size * 2);
      size_t idx = 0;
      try {
        for (; idx < prev_size; idx++) {
          new_elems[idx + prev_size] = elems_[idx];
        }
      } catch (...) {
        throw;
      }
      for (size_t j = 0; j < prev_size; j++) {
        new_elems[j] = alloc::allocate(alloc_, kBlockSize);
      }
      first_pos_ = kBlockSize * prev_size - 1;
      try {
        alloc::construct(
            alloc_,
            &new_elems[first_pos_ / kBlockSize][first_pos_ % kBlockSize],
            value);
      } catch (...) {
        first_pos_ = 0;
        throw;
      }
      last_pos_ += kBlockSize * prev_size;
      std::swap(elems_, new_elems);
    }
  }

  void push_front(T&& value) {
    if (first_pos_ >= 1) {
      first_pos_--;
      try {
        alloc::construct(
            alloc_, &elems_[first_pos_ / kBlockSize][first_pos_ % kBlockSize],
            std::move(value));
      } catch (...) {
        throw;
      }
    } else {
      size_t prev_size = elems_.size();
      std::vector<T*> new_elems(prev_size * 2);
      size_t idx = 0;
      try {
        for (; idx < prev_size; idx++) {
          new_elems[idx + prev_size] = elems_[idx];
        }
      } catch (...) {
        throw;
      }
      for (size_t j = 0; j < prev_size; j++) {
        new_elems[j] = alloc::allocate(alloc_, kBlockSize);
      }
      first_pos_ = kBlockSize * prev_size - 1;
      try {
        alloc::construct(
            alloc_,
            &new_elems[first_pos_ / kBlockSize][first_pos_ % kBlockSize],
            std::move(value));
      } catch (...) {
        first_pos_ = 0;
        throw;
      }
      last_pos_ += kBlockSize * prev_size;
      std::swap(elems_, new_elems);
    }
  }

  void pop_back() {
    if (!this->empty()) {
      alloc::destroy(
          alloc_,
          &elems_[(last_pos_ - 1) / kBlockSize][(last_pos_ - 1) % kBlockSize]);
      last_pos_--;
    }
  }

  void pop_front() {
    if (!this->empty()) {
      alloc::destroy(alloc_,
                     &elems_[first_pos_ / kBlockSize][first_pos_ % kBlockSize]);
      first_pos_++;
    }
  }

  template <class... Args>
  void emplace_back(Args&&... args) {
    if (last_pos_ == elems_.size()) {
      elems_.resize(elems_.size() * 2);
    }
    try {
      alloc::construct(alloc_,
                       &elems_[last_pos_ / kBlockSize][last_pos_ % kBlockSize],
                       std::forward<Args>(args)...);
    } catch (...) {
      throw;
    }
    ++last_pos_;
  }

  template <class... Args>
  void emplace_front(Args&&... args) {
    if (first_pos_ == 0) {
      size_t prev_size = elems_.size();
      std::vector<T*> new_elems(prev_size * 2);
      size_t idx = 0;
      try {
        for (; idx < prev_size; idx++) {
          new_elems[idx + prev_size] = elems_[idx];
        }
      } catch (...) {
        throw;
      }
      for (size_t j = 0; j < prev_size; j++) {
        new_elems[j] = alloc::allocate(alloc_, kBlockSize);
      }
      first_pos_ = kBlockSize * prev_size - 1;
      try {
        alloc::construct(
            alloc_,
            &new_elems[first_pos_ / kBlockSize][first_pos_ % kBlockSize],
            std::forward<Args>(args)...);
      } catch (...) {
        first_pos_ = 0;
        throw;
      }
      last_pos_ += kBlockSize * prev_size;
      std::swap(elems_, new_elems);
    } else {
      try {
        alloc::construct(
            alloc_, &elems_[first_pos_ / kBlockSize][first_pos_ % kBlockSize],
            std::forward<Args>(args)...);
      } catch (...) {
        throw;
      }
      --first_pos_;
    }
  }

  using iterator = common_iterator<false>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_iterator = common_iterator<true>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() { return iterator(this, 0); }

  iterator end() { return iterator(this, this->size()); }

  const_iterator cbegin() const {
    const Deque* cpy = this;
    return const_iterator(cpy, 0);
  }

  const_iterator cend() const {
    const Deque* cpy = this;
    return const_iterator(cpy, this->size());
  }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }

  reverse_iterator rbegin() { return std::reverse_iterator<iterator>(end()); }

  reverse_iterator rend() { return std::reverse_iterator<iterator>(begin()); }

  const_reverse_iterator rbegin() const {
    return std::reverse_iterator<iterator>(cend());
  }

  const_reverse_iterator rend() const {
    return std::reverse_iterator<iterator>(cbegin());
  }

  const_reverse_iterator crbegin() const {
    return std::reverse_iterator<iterator>(cend());
  }

  const_reverse_iterator crend() const {
    return std::reverse_iterator<iterator>(cbegin());
  }

  template <bool IsConst>
  void insert(common_iterator<IsConst> iter, const T& value) {
    try {
      push_back(value);
    } catch (...) {
      pop_back();
      throw;
    }
    for (auto i = --end(); i != iter; i--) {
      std::swap(*i, *(i - 1));
    }
  }

  template <bool IsConst>
  void erase(common_iterator<IsConst> iter) {
    size_t position = static_cast<size_t>(iter.pos()) + first_pos_;
    for (size_t i = position; i < last_pos_ - 1; i++) {
      elems_[i / kBlockSize][i % kBlockSize] =
          elems_[(i + 1) / kBlockSize][(i + 1) % kBlockSize];
    }
    last_pos_--;
  }

  iterator emplace(iterator position, T&& value) {
    size_t pos = position.pos() + first_pos_;
    if (pos == first_pos_ || empty()) {
      push_front(std::forward<T>(value));
      return begin();
    }
    if (pos == last_pos_) {
      push_back(std::forward<T>(value));
      return end()--;
    }
    if (pos - first_pos_ < last_pos_ - pos) {
      push_front(std::forward<T>(value));
      for (size_t i = first_pos_; i < pos; ++i) {
        elems_[i / kBlockSize][i % kBlockSize] =
            std::move(elems_[(i + 1) / kBlockSize][(i + 1) % kBlockSize]);
      }
      elems_[pos / kBlockSize][pos % kBlockSize] = std::forward<T>(value);
      return position;
    }
    push_back(std::forward<T>(value));
    for (size_t i = last_pos_ - 1; i >= pos; --i) {
      elems_[i / kBlockSize][i % kBlockSize] =
          std::move(elems_[(i - 1) / kBlockSize][(i - 1) % kBlockSize]);
    }
    elems_[pos / kBlockSize][pos % kBlockSize] = std::forward<T>(value);
    return position;
  }

  template <bool IsConst>
  class common_iterator
      : public std::iterator<std::random_access_iterator_tag, T> {
   private:
    const Deque* deque_ = nullptr;
    size_t pos_{};

   public:
    using value_type = typename std::conditional_t<IsConst, const T, T>;
    using pointer = value_type*;
    using reference = value_type&;
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = std::ptrdiff_t;

    common_iterator() = default;

    common_iterator(const Deque* deque, size_t pos)
        : pos_(pos), deque_(deque) {}

    size_t pos() const { return pos_; }

    common_iterator& operator+=(size_t n) {
      pos_ += n;
      return *this;
    }
    common_iterator& operator-=(size_t n) {
      pos_ -= n;
      return *this;
    }

    friend common_iterator operator+(const common_iterator& iter, size_t num) {
      Deque::common_iterator<IsConst> res(iter.deque_, iter.pos_);
      return (res += num);
    }
    friend common_iterator operator-(const common_iterator& iter, size_t num) {
      Deque::common_iterator<IsConst> res(iter.deque_, iter.pos_);
      return (res -= num);
    }
    friend difference_type operator-(const common_iterator& first,
                                     const common_iterator& second) {
      return (difference_type)(first.pos_ - second.pos_);
    }

    reference operator*() const {
      return deque_->elems_[(pos_ + deque_->first_pos_) / kBlockSize]
                           [(pos_ + deque_->first_pos_) % kBlockSize];
    }

    pointer operator->() const {
      return &(deque_->elems_[(pos_ + deque_->first_pos_) / kBlockSize]
                             [(pos_ + deque_->first_pos_) % kBlockSize]);
    }

    common_iterator& operator++() {
      ++pos_;
      return *this;
    }
    common_iterator operator++(int) {
      common_iterator tmp = *this;
      ++pos_;
      return tmp;
    }
    common_iterator& operator--() {
      --pos_;
      return *this;
    }
    common_iterator operator--(int) {
      common_iterator old = *this;
      --pos_;
      return old;
    }
    bool operator==(const common_iterator<IsConst>& other) const {
      return pos_ == other.pos_;
    }
    bool operator!=(const common_iterator<IsConst>& other) const {
      return pos_ != other.pos_;
    }
    bool operator<(const common_iterator<IsConst>& other) const {
      return pos_ < other.pos_;
    }
    bool operator<=(const common_iterator<IsConst>& other) const {
      return pos_ <= other.pos_;
    }
    bool operator>(const common_iterator<IsConst>& other) const {
      return pos_ > other.pos_;
    }
    bool operator>=(const common_iterator<IsConst>& other) const {
      return pos_ >= other.pos_;
    }
  };
};

