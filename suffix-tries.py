"""
Student ID: 31190863

This module consists of classes written for FIT2004 Assignment 3 S1/2021.
The classes are SDNode and SequenceDatabase for Task 1, and OrfNode and OrfFinder
for Task 2.
"""

# Task 1
class SDNode:
    """
    Node data structure for class SequenceDatabase.
    """
    def __init__(self):
        """
        Create a node for a SequenceDatabase instance. Values stored in a node is 
        the frequency of the most common sequence in `frequency`,
        `link` is a list containing the links to nodes connected 
        to the current node. Each position in `link` corresponds to a letter with 
        0 indicating the extra terminal character $ for every string:
            - A = 1, B = 2, C = 3 and D = 4
        Characters used are only [A-D].
        `data` holds the gene sequence.
        `next` holds the index of the next node that is along the path of the most 
        frequent sequence. For example, if the most frequent sequence is "ABC" then 
        the root node will hold 1, node "A" will hold 2 and so on.
        `leaf` holds the reference to the leaf node that can be reached from the current
        node that contains the most frequent sequence lexicographically
        :time complexity: best & worst: O(1) because all operations run in constant 
                          time
        :aux space complexity: O(1) because link is of constant size and
                               frequency and data contains single values
        """
        # four letters and one terminal
        self.link = [None] * 5 
        # data payload
        self.frequency = 0
        self.data = None
        self.next = None
        self.leaf = None

class SequenceDatabase:
    """
    Trie data structure implementing a sequence database for bacteria DNA sequences.
    Function addSequence() adds a sequence into the trie and query() retrieves the
    most frequent sequence in the trie.
    """
    def __init__(self):
        """
        Creates a node acting as the root to the database trie.
        :time complexity: O(1) because initialising an SDNode instance runs in 
                          constant time
        :aux space complexity: O(1) because an SDNode instance takes O(1) space
        """
        self.root = SDNode()

    def addSequence(self, s):
        """
        Adds a sequence into the sequence database using recursion. This
        function calls an auxiliary function addSequence_aux(). After adding
        s into the trie, the root node is updated if the frequency of s is 
        higher or; if root and child have equal frequencies, compare their 
        `next` value and if child is lexicographically smaller, then update 
        parent accordingly. This function does not return anything.
        :input: `s` is the new DNA sequence to be added into the database
        :time complexity: best & worst: O(len(s)) because the auxiliary function runs 
                          len(s) times
                          Calculation: O(len(s)) + O(1) = O(len(s))
        :aux space complexity: O(len(s)) because the auxiliary function uses
                               O(len(s)) space. compare_freq() takes O(1) space
                               Calculation: O(len(s)) + O(1) = O(len(s))
        """
        current = self.root
        # add the sequence into the trie
        self.addSequence_aux(current, s)              # runs in O(len(s)) time
        index = ord(s[0]) - 65 + 1
        # get the node for the first character in s
        child_node = current.link[index]
        # update the root data if necessary
        self.compare_freq(child_node, index, current) # runs in O(1) time
    
    def addSequence_aux(self, current, s, counter=0):
        """
        This function goes through each character in the s and updates the frequency 
        of the leaf node if s exists. If not, a path to a new leaf node is created 
        and the frequency is set to 1. When the recursion returns, it returns a 
        reference of the leaf node and the index of the next character the current node 
        points to. If the parent node has a frequency less than the child node or if the 
        frequency is the same and the child is lexicographically smaller, then the leaf
        node reference, `next` value and frequency of the parent node will be recursively 
        updated. 
        To compare lexicographically, only the `next` value of the two child nodes on the
        same level are compared. This happens only when two nodes are children of one parent
        node. This is because when traversing from the root up till the parent, the 
        sequence we are adding will have the same prefix as the current most frequent 
        sequence up till a node with more than one child.
        At the end of execution of this function, each node will contain `next` value,
        the frequency and data (sequence) of the most frequent and lexicographically 
        small sequence that has a prefix from the root up till that node.
        :input: `current` is the root node of the database;
                `s` is the sequence to be added into the database;
                `counter` keeps track of the number of characters already processed
        :time complexity: best & worst: O(len(s)) because s needs to be iterated through 
                          to the leaf and back to the first character
                          Calculation: O(len(s)) + O(len(s)) + O(len(s)) + O(1) = O(len(s))
                          Note: Creation of nodes are O(1) and string comparison is 
                                done by integer comparison using `next` values. Copying
                                data to the leaf node takes O(len(s)) time
        :aux space complexity: O(len(s)) as the depth of the recursive stack is of
                               len(s) and data stored in the leaf is O(len(s)). Nodes 
                               store single values so it is O(1). compare_freq() uses O(1) 
                               space
                               Calculation: O(len(s)) + O(len(s)) + O(len(s)) + O(1) = O(len(s))
        """
        # if it the last letter in the string (terminal character)
        if len(s) == counter: 
            # if linked node (leaf) does not exist then create a linked node
            if current.link[0] is None:
                current.link[0] = SDNode()
            # update the frequency of the leaf
            current.link[0].frequency += 1
            # if the frequency of the leaf is higher, then update the parent's next value
            if current.frequency < current.link[0].frequency:
                current.next = 0
            # update new leaf node
            current = current.link[0]
            current.data = s        # runs in O(len(s)) - copying data
            current.next = 0
            current.leaf = current
            return (current, 0)
        else:
            # calculate the index of the character
            index = ord(s[counter]) - 65 + 1
            # if linked node exists then traverse to linked node
            if current.link[index] is not None:
                current = current.link[index]
            # if not, create a new linked node and traverse to it
            else:
                current.link[index] = SDNode()
                current = current.link[index]   
            # increment counter to keep count of characters processed
            counter += 1
            # recursively add characters into the trie
            child = self.addSequence_aux(current, s, counter)   # runs in O(len(s)) time
            self.compare_freq(child[0], child[1], current)      # runs in O(1) time
            # index is the `next` value for the parent node if an update is required
            return (current, index)
    
    def compare_freq(self, child, nxt, parent):
        """
        This function compares the frequency of the child and parent node. If the child
        node has a higher frequency then the parent's data payload will be updated.
        If parent and child have the same frequency, then update parent if the child is
        lexicographically smaller
        :input: child node, nxt (an integer representing the index of the current node character),
                parent node
        :time complexity: best & worst: O(1) because all assigning of values is done in constant time
        :aux space complexity: O(1) because all variables store single values
        """
        # compare frequency of child node and parent node
        # and update parent if child node has higher frequency
        if child.frequency > parent.frequency:
            parent.leaf = child.leaf
            parent.frequency = child.frequency
            parent.next = nxt
        # if parent and child have equal frequencies, compare their `next` value
        # if there is a None or the child is lexicographically smaller, then update parent
        elif child.frequency == parent.frequency:
            if nxt is None or nxt <= parent.next:
                parent.leaf = child.leaf
                parent.frequency = child.frequency
                parent.next = nxt

    def query(self, q):
        """
        This function searches the database and returns the most frequent and lexicographically
        small sequence in the database with the prefix `q`. If q is empty then return the data
        stored at the root of the trie. If there are no strings in the database, the function
        will return None.
        :input: a prefix of a sequence to be searched for in the database, q
        :output: the most common and lexicographically small sequence in the database with prefix
                 q; if the prefix does not exist in the database, function will return None
        :time complexity: best & worst: O(len(q)) because function needs to go loop through every 
                          character in `q` and traverse the database to the node of that character
        :aux space complexity: O(1) because `current` only saves the current node
        """
        current = self.root
        if len(q) == 0:     # if the query string is empty
            if current.leaf is None:
                return current.leaf        # return None if the database is empty
            else:
                return current.leaf.data   # return the most frequent string in the database
        for char in q:      # runs in O(len(q)) time
            index = ord(char) - 65 + 1
            if current.link[index] is not None:
                current = current.link[index]
            else:
                return None
        return current.leaf.data

# Task 2
class OrfNode:
    """
    Node data structure for OrfFinder trie class.
    """
    def __init__(self, g=None):
        """
        Creates a node for a OrfFinder instance. Values stored include a list, 
        `position` which stores the index of the character in the original genome 
        string that has a path through the node when a using a suffix to traverse 
        the trie. For example, in AAAB, the first node linked to the root will 
        store position = [0,1,2] because suffixes AAAB, AAB, AB will pass through 
        the node.
        `link` is a list containing the links to nodes connected 
        to the current node. Each position in `link` corresponds to a letter with 
        0 indicating the extra terminal character $ for every string:
            - A = 1, B = 2, C = 3 and D = 4
        Characters used are only [A-D].
        `genome` stores the original genome string used to construct the trie. Only
        the root node will store this genome. By default, all other nodes will store
        None in this attribute.
        :input: g - default value is None; for the root, the genome string will be 
                passed in
        :time complexity: best O(1) when `genome` is None
                          worst O(len(g)) because g needs to be copied into `genome`
        :aux space complexity: O(len(g)) if `genome` stores the string g
                               O(1) because link array is constant size and position
                               and genome store single values
        """
        self.link = [None]*5
        # data payload
        self.position = []
        self.genome = g         # assign the default value None

class OrfFinder:
    """
    OrfFinder class (constructed as a suffix trie) to store a genome string.
    Function find() is used to find all substrings in the genome string using the trie
    that start and end with certain characters.
    """

    def __init__(self, genome):
        """
        Creates a suffix trie for the genome string of length N+1 including a terminal symbol. 
        The terminal symbol is lexicographically smaller than A.
        :input: genome - a string consisting of uppercase letters [A-D] of length N
        :time complexity: best & worst O(N^2) because the nested loop to add suffixes to the 
                          trie runs in O(N^2) time. 
                          Creating nodes and accessing the linked nodes is done in constant time.
                          Creating the root node is O(N) time
                          Calculation: O(N) + O(N^2) = O(N^2)
        :aux space complexity: O(N^2) as we store all suffixes in the trie
        """
        # create root and store the genome string
        self.root = OrfNode(genome)     # runs in O(N) time
        j = 0
        # the loop below runs in O(N^2) time
        while j < len(genome):          # runs in O(N) time
            current = self.root         # reset current to the root before adding each suffix
            for i in range(j, len(genome)+1):       # runs in O(N) time
                if i == len(genome):    # at the terminal symbol (end of the string)
                    if current.link[0] is not None:
                        current = current.link[0]
                    else:       # create a node if link has no node
                        current.link[0] = OrfNode()             # runs in O(1) time
                        current = current.link[0]
                else:
                    index = ord(genome[i]) - 65 + 1
                    if current.link[index] is not None:
                        current = current.link[index]
                    else:       # create a node if link has no node
                        current.link[index] = OrfNode()         # runs in O(1) time
                        current = current.link[index]
                # store the position of the current character in the node data payload
                current.position.append(i)
            j += 1
    
    def find(self, start, end):
        """
        Finds all substrings in the original genome string that starts with `start` and
        ends with `end`. `start` and `end` do not overlap. Both are non-empty strings.
        :input: start - a string that is to be the start of each substring;
                end - a string that is to be the end of each substring
        :time complexity: best & worst: O(len(start) +  U + len(end))
                          U is the number of characters in the output list
        :aux space complexity: O(U) because the return list stores O(U) characters
        """
        current = self.root
        # loop through start to ensure that it exists in the genome string
        for char in start:      # runs in O(len(start))
            index = ord(char) - 65 + 1
            if current.link[index] is not None:
                current = current.link[index]
            else:               # if start does not exist in the genome string
                return []
        # indexes refer to the last index (end) of the `start` string in the genome string
        # if start = "AA" and genome = "AAB", start_index = [1]
        start_index = current.position      # list of indexes
        # traverse from the root again to find end since it is a suffix trie
        current = self.root 
        for char in end:
            index = ord(char) - 65 + 1
            if current.link[index] is not None:
                current = current.link[index]
            else:               # if start does not exist in the genome string
                return []
        # indexes refer to the last index (end) of the `end` string in the genome string
        end_index = current.position
        substrings = []
        genome = self.root.genome       # get genome string stored in root

        # the loop below runs in U time because of the list slicing;

        # len(start_index)+len(end_index) is less than U if either start or
        # end has length > 1 so the loop will run in < U time but the slicing takes
        # O(U) time so the block of code will run in O(U) time

        # if len(start) and len(end) is 1,
        # len(start_index) and len(end_index) are both at most of length len(genome) each
        # and are of len(genome) each if the genome is a uniform string
        # In this case, the length of the substring is at least 2 so the slicing will 
        # run in N^2 time or U time which is equal to the complexity of the nested 
        # loop which makes it O(N^2 + N^2) = O(N^2) and since U = N^2, this loop
        # below runs in U time
        for index in start_index:  
            for idx in end_index:
                last_char = idx - (len(end) - 1)
                if idx <= index or index >= last_char:  # start and end cannot overlap
                    continue
                else:
                    # since start_index contains the last index, starting index is calculated
                    first_char = index - (len(start)-1)
                    substrings.append(genome[first_char:idx+1])
        return substrings


