### Документация для класса Deque

#### Описание
`Deque` (двусторонняя очередь) представляет собой контейнер, который обеспечивает эффективное добавление и удаление элементов как в начале, так и в конце контейнера.

### Методы

1. **push_back(value):**
   - Добавляет элемент с заданным значением в конец деки.

```cpp
void push_back(const T& value);
```

2. **push_front(value):**
   - Добавляет элемент с заданным значением в начало деки.

```cpp
void push_front(const T& value);
```

3. **pop_back():**
   - Удаляет последний элемент из деки.

```cpp
void pop_back();
```

4. **pop_front():**
   - Удаляет первый элемент из деки.

```cpp
void pop_front();
```

5. **at(index):**
   - Возвращает ссылку на элемент по указанному индексу в деке.

```cpp
T& at(size_t index);
```

6. **insert(iter, value):**
   - Вставляет элемент со значением `value` перед элементом, на который указывает итератор `iter`.

```cpp
void insert(iterator iter, const T& value);
```

7. **erase(iter):**
   - Удаляет элемент, на который указывает итератор `iter`.

```cpp
void erase(iterator iter);
```

8. **size():**
   - Возвращает количество элементов в деке.

```cpp
size_t size() const;
```

9. **empty():**
   - Проверяет, пуста ли дека.

```cpp
bool empty() const;
```

### Конструкторы

1. **Deque():**
   - Конструктор по умолчанию.

```cpp
Deque() = default;
```

2. **Deque(copy):**
   - Конструктор копирования.

```cpp
Deque(const Deque<T>& copy);
```

3. **Deque(count):**
   - Конструктор, создающий деку с `count` элементами.

```cpp
Deque(const size_t& count);
```

4. **Deque(count, value):**
   - Конструктор, создающий деку с `count` элементами, каждый из которых равен `value`.

```cpp
Deque(const size_t& count, const T& value);
```

### Деструктор

```cpp
~Deque();
```

### Операторы

1. **operator=():**
   - Оператор присваивания.

```cpp
void operator=(const Deque<T>& copy);
```

2. **operator[]():**
   - Возвращает ссылку на элемент по указанному индексу в деке.

```cpp
T& operator[](const size_t& index);
```

### Итераторы

1. **begin():**
   - Возвращает итератор, указывающий на первый элемент деки.

```cpp
iterator begin();
```

2. **end():**
   - Возвращает итератор, указывающий на элемент, следующий за последним элементом деки.

```cpp
iterator end();
```

3. **cbegin():**
   - Возвращает константный итератор, указывающий на первый элемент деки.

```cpp
const_iterator cbegin();
```

4. **cend():**
   - Возвращает константный итератор, указывающий на элемент, следующий за последним элементом деки.

```cpp
const_iterator cend();
```

5. **rbegin():**
   - Возвращает обратный итератор, указывающий на последний элемент деки.

```cpp
reverse_iterator rbegin() noexcept;
```

6. **rend():**
   - Возвращает обратный итератор, указывающий на элемент перед первым элементом деки.

```cpp
reverse_iterator rend() noexcept;
```

7. **crbegin():**
   - Возвращает константный обратный итератор, указывающий на последний элемент деки.

```cpp
const_reverse_iterator crbegin() const noexcept;
```

8. **crend():**
   - Возвращает константный обратный итератор, указывающий на элемент перед первым элементом деки.

```cpp
const_reverse_iterator crend() const noexcept;
```

### Утилиты

1. **get_elem_count():**
   - Возвращает количество элементов в деке.

```cpp
size_t get_elem_count();
```

2. **cout_vector():**
   - Выводит содержимое деки в консоль в виде блоков.

```cpp
void cout_vector();
```
