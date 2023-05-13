#include <algorithm>
#include <iostream>
#include <vector>

template <typename T>
class Deque {
 private:
  const size_t kStartingSize = 2;  // 2
  const size_t kArraySize = 10;    // 10
  std::vector<T*> blocks_;
  int64_t start_ = 0;
  size_t finish_ = 0;
  size_t starting_block_ = 0;
  size_t ending_block_ = 0;
  size_t size_ = 0;

 public:
  // Constructors
  Deque() = default;

  // might need a change
  Deque(const Deque<T>& copy)
      : size_(copy.size_),
        start_(copy.start_),
        finish_(copy.finish_),
        starting_block_(copy.starting_block_),
        ending_block_(copy.ending_block_) {
    blocks_.resize(copy.blocks_.size());
    for (size_t i = 0; i < copy.blocks_.size() * kArraySize; ++i) {
      if (blocks_[i / kArraySize] == nullptr) {
        blocks_[i / kArraySize] =
            reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize));
      }
      blocks_[i / kArraySize][i % kArraySize] =
          copy.blocks_[i / kArraySize][i % kArraySize];
    }
  }

  Deque(const size_t& count) : size_(count) {
    size_t capacity =
        std::max(kStartingSize, (kStartingSize * count) / kArraySize + 1);
    blocks_.resize(capacity);

    size_t el_count = count;
    for (size_t i = 0; i < capacity; ++i) {
      blocks_[i] = reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize));
      if constexpr (std::is_default_constructible_v<T> &&
                    std::is_constructible_v<T>) {
        for (size_t j = 0; j < kArraySize; ++j) {
          if (el_count-- != 0) {
            new (&blocks_[i][j]) T();
          }
        }
      } else {
      }
    }

    starting_block_ = 0;
    start_ = 0;
    ending_block_ = count / kArraySize;
    finish_ = count % kArraySize;
  }

  Deque(const size_t& count, const T& value) : size_(count) {
    starting_block_ = 0;
    start_ = 0;
    ending_block_ = count / kArraySize;
    finish_ = count % kArraySize;
    try {
      size_t capacity =
          std::max(kStartingSize, (kStartingSize * count) / kArraySize + 1);
      blocks_.resize(capacity);

      for (size_t i = 0; i < capacity; ++i) {
        blocks_[i] =
            reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize));
      }

      capacity *= kArraySize;

      size_t el_count = count;
      for (size_t i = 0; i < size_; ++i) {
        if (el_count-- != 0) {
          new (&blocks_[i / kArraySize][i % kArraySize]) T(value);
        }
      }
    } catch (...) {
      for (size_t i = 0; i < blocks_.size(); ++i) {
        operator delete[](blocks_[i]);
      }
      throw 1;
    }
  }

  // Destructor
  ~Deque() {
    for (size_t i = 0; i < blocks_.size(); ++i) {
      operator delete[](blocks_[i]);
    }
  }

  // Operators
  void operator=(const Deque<T>& copy) {
    Deque<T> temp(copy);
    std::swap(blocks_, temp.blocks_);
    std::swap(starting_block_, temp.starting_block_);
    std::swap(ending_block_, temp.ending_block_);
    std::swap(start_, temp.start_);
    std::swap(finish_, temp.finish_);
    std::swap(size_, temp.size_);
  }

  T& operator[](const size_t& index) {
    size_t actual_index = index + kArraySize * starting_block_ + start_;
    return blocks_[actual_index / kArraySize][actual_index % kArraySize];
  }

  const T& operator[](const size_t& index) const {
    size_t actual_index = index + kArraySize * starting_block_ + start_;
    return blocks_[actual_index / kArraySize][actual_index % kArraySize];
  }

  // Getters
  size_t size() const { return size_; }

  bool empty() const { return size_ == 0; }

  // Methods
  void push_back(const T& value) {
    if (empty()) {
      blocks_.resize(
          1, reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize)));
      new (&blocks_[0][0]) T(value);
      ++size_;
      ++finish_;
    } else if (finish_ < kArraySize) {
      new (&blocks_[ending_block_][finish_]) T(value);
      ++size_;

      ++finish_;
    } else if (ending_block_ < blocks_.size() - 1) {
      finish_ = 0;
      ++ending_block_;
      new (&blocks_[ending_block_][finish_]) T(value);
      ++size_;
      ++finish_;
    } else {
      size_t initial_size = blocks_.size();
      blocks_.resize(initial_size * 3);
      for (size_t i = 0; i < initial_size; ++i) {
        blocks_[i + initial_size] = blocks_[i];
        blocks_[i + 2 * initial_size] =
            reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize));
        blocks_[i] =
            reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize));
      }
      starting_block_ += initial_size;
      ending_block_ += initial_size;
      this->push_back(value);
    }
  }

  void push_front(const T& value) {
    if (empty()) {
      blocks_.resize(
          1, reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize)));
      new (&blocks_[0][0]) T(value);
      ++size_;
      ++finish_;
    } else if (start_ - 1 >= 0) {
      --start_;
      ++size_;
      new (&blocks_[starting_block_][start_]) T(value);
    } else if (starting_block_ > 0) {
      start_ = kArraySize - 1;
      --starting_block_;
      ++size_;
      new (&blocks_[starting_block_][start_]) T(value);
    } else {
      size_t initial_size = blocks_.size();
      blocks_.resize(initial_size * 3);
      for (size_t i = 0; i < initial_size; ++i) {
        blocks_[i + initial_size] = blocks_[i];
        blocks_[i + 2 * initial_size] =
            reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize));
        blocks_[i] =
            reinterpret_cast<T*>(operator new[](sizeof(T) * kArraySize));
      }
      starting_block_ += initial_size;
      ending_block_ += initial_size;
      this->push_front(value);
    }
  }

  void pop_back() {
    if (!empty()) {
      --finish_;
      blocks_[ending_block_][finish_].~T();
      if (finish_ == 0) {
        finish_ = kArraySize;
        --ending_block_;
      }
      --size_;
    }
  }

  void pop_front() {
    if (!empty()) {
      blocks_[starting_block_][start_].~T();
      ++start_;
      if (start_ == int64_t(kArraySize)) {
        start_ = 0;
        ++starting_block_;
      }
      --size_;
    }
  }

  T& at(size_t index) {
    if (index >= get_elem_count()) {
      throw std::out_of_range("out_of_range");
    }
    size_t actual_index = index + kArraySize * starting_block_ + start_;
    return blocks_[actual_index / kArraySize][actual_index % kArraySize];
  }

  const T& at(size_t index) const {
    if (index >= get_elem_count()) {
      throw std::out_of_range("out of range");
    }
    size_t actual_index = index + kArraySize * starting_block_ + start_;
    return blocks_[actual_index / kArraySize][actual_index % kArraySize];
  }

  // Utility
  size_t get_elem_count() {
    return ending_block_ * kArraySize + finish_ - starting_block_ * kArraySize -
           start_;
  }

  void cout_vector() {
    if (get_elem_count() != size_) {
      std::cout << "SIZES DO NOT MATCH" << std::endl;
      std::cout << ending_block_ * kArraySize + finish_ -
                       starting_block_ * kArraySize - start_
                << "!=" << size_ << std::endl;
    }
    for (size_t i = 0; i < blocks_.size(); ++i) {
      std::cout << '[' << i << ']' << " : ";
      if (blocks_[i] != nullptr) {
        for (size_t j = 0; j < kArraySize; ++j) {
          if (blocks_[i][j] >= 0) {
            std::cout << blocks_[i][j] << " ";
          } else {
            std::cout << '#' << " ";
          }
        }
        std::cout << std::endl;
      }
    }
    std::cout << "----------------------------------" << std::endl;
  }

  // Iterators
  template <bool IsConst = false>
  class Iterator
      : std::iterator<std::random_access_iterator_tag,
                      typename std::conditional<IsConst, const T, T>::type> {
   private:
    Deque<T>* pointer_to_deque_;
    size_t index_ = 0;

   public:
    using value_type = typename std::iterator<
        std::random_access_iterator_tag,
        typename std::conditional<IsConst, const T, T>::type>::value_type;
    using difference_type = typename std::iterator<
        std::random_access_iterator_tag,
        typename std::conditional<IsConst, const T, T>::type>::difference_type;
    using iterator_category =
        typename std::iterator<std::random_access_iterator_tag,
                               typename std::conditional<IsConst, const T, T>::
                                   type>::iterator_category;
    using pointer = typename std::iterator<
        std::random_access_iterator_tag,
        typename std::conditional<IsConst, const T, T>::type>::pointer;
    using reference = typename std::iterator<
        std::random_access_iterator_tag,
        typename std::conditional<IsConst, const T, T>::type>::reference;

    // Constructors
    Iterator() = default;

    Iterator(Deque<T>* deque, size_t index)
        : pointer_to_deque_(deque), index_(index) {}

    Iterator(const Iterator<IsConst>& copy)
        : pointer_to_deque_(copy.pointer_to_deque_), index_(copy.index_) {}

    // Utility
    Deque<T>* get_deque() { return pointer_to_deque_; }

    size_t get_index() { return index_; }

    // Operators
    void operator=(const Iterator& copy) {
      index_ = copy.index_;
      pointer_to_deque_ = copy.pointer_to_deque_;
    }
    Iterator<IsConst>& operator++() {
      ++index_;
      return *this;
    }

    Iterator<IsConst> operator++(int) {
      Iterator<IsConst> temp_iter(*this);
      ++index_;
      return temp_iter;
    }

    Iterator<IsConst>& operator--() {
      --index_;
      return *this;
    }

    Iterator<IsConst> operator--(int) {
      Iterator<IsConst> temp_iter(*this);
      --index_;
      return temp_iter;
    }

    Iterator<IsConst>& operator+=(difference_type iter_jump) {
      index_ += iter_jump;
      return *this;
    }

    Iterator<IsConst>& operator-=(difference_type iter_jump) {
      index_ -= iter_jump;
      return *this;
    }

    Iterator<IsConst> operator+(difference_type iter_jump) const {
      Iterator<IsConst> temp_iter(*this);
      temp_iter.index_ += iter_jump;
      return temp_iter;
    }

    Iterator<IsConst> operator-(difference_type iter_jump) const {
      Iterator<IsConst> temp_iter(*this);
      temp_iter.index_ -= iter_jump;
      return temp_iter;
    }

    difference_type operator-(const Iterator<IsConst>& iter) const {
      return index_ - iter.index_;
    }

    reference operator*() const { return (*pointer_to_deque_)[index_]; }

    pointer operator->() const { return &((*pointer_to_deque_)[index_]); }

    friend bool operator==(const Iterator<IsConst>& first,
                           const Iterator<IsConst>& second) {
      return (first.pointer_to_deque_ == second.pointer_to_deque_ &&
              first.index_ == second.index_);
    }

    friend bool operator>(const Iterator<IsConst>& first,
                          const Iterator<IsConst>& second) {
      return (first.pointer_to_deque_ == second.pointer_to_deque_ &&
              first.index_ > second.index_);
    }

    friend bool operator!=(const Iterator<IsConst>& first,
                           const Iterator<IsConst>& second) {
      return !(first == second);
    }

    friend bool operator<(const Iterator<IsConst>& first,
                          const Iterator<IsConst>& second) {
      return (first != second && !(first > second));
    }

    friend bool operator<=(const Iterator<IsConst>& first,
                           const Iterator<IsConst>& second) {
      return !(first > second);
    }

    friend bool operator>=(const Iterator<IsConst>& first,
                           const Iterator<IsConst>& second) {
      return !(first < second);
    }
  };

  // Methods to work with iterators
  using iterator = Iterator<false>;
  using const_iterator = Iterator<true>;
  using reverse_iterator = std::reverse_iterator<Iterator<false>>;
  using const_reverse_iterator = std::reverse_iterator<Iterator<false>>;

  iterator begin() { return iterator(this, 0); }

  iterator end() { return iterator(this, size_); }

  const_iterator cbegin() { return const_iterator(this, 0); }

  const_iterator cend() { return const_iterator(this, size_); }

  reverse_iterator rbegin() noexcept {
    return reverse_iterator(iterator(this, size_));
  }

  reverse_iterator rend() noexcept {
    return reverse_iterator(iterator(this, 0));
  }

  const_reverse_iterator crbegin() const noexcept {
    return const_reverse_iterator(iterator(this, size_));
  }

  const_reverse_iterator crend() const noexcept {
    return const_reverse_iterator(iterator(this, 0));
  }

  void insert(iterator iter, const T& value) {
    push_back(value);
    for (size_t i = size() - 1; i > iter.get_index(); --i) {
      std::swap(this->at(i), this->at(i - 1));
    }
  }

  void erase(iterator iter) {
    for (size_t i = iter.get_index(); i < size() - 1; ++i) {
      std::swap(this->at(i), this->at(i + 1));
    }
    if (!empty()) {
      pop_back();
    }
  }
};
