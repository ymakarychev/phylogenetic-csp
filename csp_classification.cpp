#include <bitset>
#include <string>
#include <array>
#include <algorithm>
#include <ranges>
#include <iostream>

using Permutation = std::array<int, 4>;

enum class TreeType {
    C, P
};

enum class TreeName {
    C01, C02, C03, C10, C12, C13, C20, C21, C23, C30, C31, C32, P1, P2, P3
};

constexpr int tree_count = 15;
constexpr int constraint_count = 1 << tree_count;

std::pair<int, int> other_two(int a, int b) {
    int c = -1;
    int d = -1;

    for (int i = 0; i < 4; ++i) {
        if (i != a && i != b) {
            if (c == -1) c = i;
            else d = i;
        }
    }

    return {c, d};
}

class Tree {
    TreeName id;    
    
    public:
    constexpr Tree(TreeName id) : id(id) {};
    constexpr Tree(int id) : id(static_cast<TreeName>(id)) {};
    
    constexpr TreeName get_tree_name() const {
        return id;
    }
    constexpr int get_id() const {
        return static_cast<int>(id);
    }

    constexpr TreeType get_type() const {
        return (get_id() < static_cast<int>(TreeName::P1)) ? TreeType::C : TreeType::P;
    }
    
    static constexpr Tree get_C_tree(int i, int j) {        
        return Tree(i*3 + j - (i < j ? 1 : 0));
    }

    static constexpr Tree get_P_tree(int i) {
        return Tree(static_cast<TreeName>(static_cast<int>(TreeName::P1) - 1 + i));
    }

    // Use indices indexed from 1 to 4 in the output
    inline static const std::array<std::string, 15> tree_names = {
        "C12", "C13", "C14", "C21", "C23", "C24", "C31", "C32", "C34", "C41", "C42", "C43",
        "P2", "P3", "P4"
    };

    std::string to_string() const {
        return tree_names[static_cast<int>(id)];
    }

    constexpr Tree apply_permutation(const Permutation& perm) const {
        if (get_type() == TreeType::C) {
            int i = get_id() / 3;
            int j = get_id() % 3;
            if (j >= i) j++;
            
            return get_C_tree(perm[i], perm[j]);
        } else {
            int i = get_id() - static_cast<int>(TreeName::P1) + 1;
            if (perm[0] == 0) {
                return get_P_tree(perm[i]);
            }

            if (perm[i] == 0) {
                return get_P_tree(perm[0]);
            }
            
            return get_P_tree(6 - perm[0] - perm[i]); // exclude perm[0] and perm[1] from the sum 1 + 2 + 3
        }
    }
};

struct permutation_action{
    std::array<Permutation, 24> all_permutations;
    std::array<std::array<int, 24>, tree_count> tree_orbits;

    permutation_action() {
        Permutation perm = {0, 1, 2, 3};
        for (int i = 0; i < 24; i++) {
            all_permutations[i] = perm;
            std::next_permutation(perm.begin(), perm.end());
        }

        for (int i = 0; i < tree_count; ++i) {
            Tree tree(i);
            for (int j = 0; j < 24; j++) {
                tree_orbits.at(i).at(j) = tree.apply_permutation(all_permutations[j]).get_id();
            }
        }
    }
};

class Constraint {
    std::bitset<tree_count> trees;
    inline static permutation_action action;

    public:       
    Constraint() : trees(0) {};
    Constraint(std::bitset<tree_count> trees) : trees(trees) {};
    Constraint(int csp_id) : trees(csp_id) {};

    int get_csp_id() const {
        return static_cast<int>(trees.to_ulong());
    }

    bool contains_tree(const Tree tree) {
        return trees[tree.get_id()];
    }

    Constraint& add(const TreeName tree) {
        trees[static_cast<int>(tree)] = 1;
        return *this;
    }

    Constraint& add(int tree_id) {
        trees[tree_id] = 1;
        return *this;
    }

    static Constraint get_all_caterpillars()
    {
        return Constraint((1 << static_cast<int>(TreeName::P1)) - 1);        
    } 

    auto get_orbit() const{
        return std::views::iota(0, 24) | std::views::transform([this](int i) {
            Constraint target;
            for (int j = 0; j < tree_count; ++j) {
                if (trees[j]) {
                    target.add(action.tree_orbits[j][i]);
                }
            }
            return target; 
        });
    }

    std::string to_string() const {
        std::string str_repr = "{";
        bool first = true;
        for (int i = 0; i < tree_count; i++) {
            if (trees[i]) {
                if (!first) str_repr += ", ";
                str_repr += Tree(i).to_string();
                first = false;
            }
        }
        str_repr += "}";
        return str_repr;
    }

    // returns whether a is a subset of b
    friend bool subset(const Constraint& a, const Constraint& b) {
        return (a.trees & b.trees) == a.trees;
    }

    friend Constraint intersection(const Constraint& a, const Constraint& b) {
        return Constraint(a.trees & b.trees);
    }

    friend Constraint operator +(Constraint &a, const Constraint& b) {        
        return Constraint(a.trees | b.trees);
    }
};

struct ConstraintInfo
{
    Constraint constraint = Constraint();
    bool poly_tractable = false;
    enum {Unknown, Canonical, Noncanonical} canonical_representative = Unknown;
    bool generator = false;
    bool always_satisfiable = false;
    int relaxation_id = -1;    
};

class ConstraintTable {    
    public:
    std::array<ConstraintInfo, constraint_count> table;
    ConstraintTable() {
        for (int i = 0; i < constraint_count; i++) {
            table[i].constraint = Constraint(i);
        }
    }

    void add_generator(const Constraint& constr) {
        table[constr.get_csp_id()].generator = true;
        table[constr.get_csp_id()].poly_tractable = true;
    }

    void add_generator_orbit(const Constraint& base_constr) {
        for (auto c : base_constr.get_orbit()) {
            add_generator(c);
        }
    }

    void set_generators_4_1100(){
        // set generators with linear space orthogonal to (1,1,0,0)
        // that is, 0 and 1 must be on the same side at the top level
        // that is neither of them can be the first in a C tree
        Constraint base_constr;
        base_constr
            .add(TreeName::C20)
            .add(TreeName::C21)
            .add(TreeName::C23)
            .add(TreeName::C30)
            .add(TreeName::C31)
            .add(TreeName::C32)
            .add(TreeName::P1);            

        add_generator_orbit(base_constr);
    }


    void set_generators_4_1111(){ 
        // set generators with linear space orthogonal to (1,1,1,1)
        // this is an arbitrary 2-2 split
        Constraint base_constr;
        base_constr
            .add(TreeName::P1)
            .add(TreeName::P2)
            .add(TreeName::P3);

        add_generator(base_constr); // the orbit is trivial, sufficient to add the base constraint only
    }

    void set_generators_3_110(){ 
        // set generators with linear space orthogonal to (1,1,0), element 3 is ignored
        // that is, the first split that separates some of 0, 1, and 2 must keep 0 and 1 together
        // if this is a C tree, then the top elements should be (2, anything) or (3, 2)
        // the only P tree is P1
        Constraint base_constr;
        base_constr                         
            .add(TreeName::C20)
            .add(TreeName::C21)
            .add(TreeName::C23)
            .add(TreeName::C32)
            .add(TreeName::P1);

        add_generator_orbit(base_constr);
    }

    int compute_poly_tractable_constraints() {
        int count = 0;          //the count does not include the generators
        set_generators_4_1100();
        set_generators_4_1111();
        set_generators_3_110();

        for (int i = constraint_count - 1; i >= 0; --i) {
            if (table[i].poly_tractable) {
                for (int j = constraint_count - 1; j > i; --j) {
                    if (table[j].poly_tractable) {
                            int r = intersection(Constraint(i), Constraint(j)).get_csp_id();
                            if (!table[r].poly_tractable) {
                                table[r].poly_tractable = true;
                                count++;                        
                        }}}}}
        return count;
    }

    void mark_canonical_representatives() {
        for (int i = 0; i < constraint_count; i++) {
            if (table[i].canonical_representative == ConstraintInfo::Unknown) {
                table[i].canonical_representative = ConstraintInfo::Canonical;
                for (auto c : table[i].constraint.get_orbit()) {
                    if (c.get_csp_id() != i) {
                        table[c.get_csp_id()].canonical_representative = ConstraintInfo::Noncanonical;
                    }
                }
            }
        }
    }

    void mark_always_satisfiable() {
        auto mark_orbit = [this](Constraint& base_constr) {
            for (auto c : base_constr.get_orbit()) {
                table[c.get_csp_id()].always_satisfiable = true;
            }
        };
        
        auto cat = Constraint::get_all_caterpillars();
        mark_orbit(cat);
        cat = cat.add(TreeName::P1);
        mark_orbit(cat);
        cat.add(TreeName::P2);
        mark_orbit(cat);
        cat.add(TreeName::P3);
        mark_orbit(cat);
    }

    void compute_relaxations() {
        for (int csp = 0; csp < constraint_count; csp++) {
            int relaxation = csp;
            while (relaxation < constraint_count) {
                
                if (table[relaxation].poly_tractable && subset(Constraint(csp), Constraint(relaxation))) {
                    table[csp].relaxation_id = relaxation;
                    break;
                }
                relaxation += 1;
            }
        }
    }   

    void process_all(){
        mark_canonical_representatives();
        compute_poly_tractable_constraints();
        compute_relaxations();
    }
};

int main() {
    std::cout << "Phylogenetic CSPs classification.\n";
        
    ConstraintTable csp_table;    
    csp_table.process_all();

    //print information about the poly-tractable classes
    int isomorphism_classes = 0;
    int poly_tractable_classes = 0;
    for (auto record : csp_table.table) {
        if (record.canonical_representative == ConstraintInfo::Canonical) {
            isomorphism_classes++;
            if (record.poly_tractable) {
                poly_tractable_classes++;
                std::cout << poly_tractable_classes << ". Canonical representative: " << record.constraint.to_string() << "\n";
            } 
        }        
    }

    //print information about the always satisfiable classes
    Constraint cat = Constraint::get_all_caterpillars();
    std::cout << "\nAdditionally, all instances of the following phylogenetic CSPs are always satisfiable and thus poly-tractable:\n" << std::flush;
    std::cout << "a. " << cat.to_string() << "\n";
    cat = cat.add(TreeName::P1);
    std::cout << "b. " << cat.to_string () << "\n";
    cat.add(TreeName::P2);
    std::cout << "c. " << cat.to_string() << "\n";
    cat.add(TreeName::P3);
    std::cout << "d. " << cat.to_string() << "\n\n";

    //print information about relaxations
    std::cout << "The following constraints have non-trivial relaxations (only one representative for each CSP class is printed):\n";
    int relaxation_count = 0;
    for (auto record : csp_table.table) {
        if (record.canonical_representative == ConstraintInfo::Canonical) {            
            if (!record.poly_tractable && record.relaxation_id != -1) {
                relaxation_count++;
                std::cout << relaxation_count << ". Relaxation for " << record.constraint.to_string() << " is " << Constraint(record.relaxation_id).to_string() 
                << " Delta: " << Constraint(record.relaxation_id - record.constraint.get_csp_id()).to_string() << "\n";
            }
        }
    }    

    //print statistics
    std::cout << "Total constraints: " << constraint_count << "\n" << std::flush;
    std::cout << "Isomorphism classes: " << isomorphism_classes << "\n";
    std::cout << "Poly-tractable classes, excluding always satisfiable: " << poly_tractable_classes << "\n";   
    std::cout << "Always satisfiable classes: " << 4 << "\n"; 
    std::cout << "Constraints with non-trivial relaxations: " << relaxation_count << "\n";
    return 0;  
}
