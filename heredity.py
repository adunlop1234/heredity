import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}

PARENT_PROB = {
    
    "parent": {

        0: {
            0: 1 - PROBS["mutation"],
            1: PROBS["mutation"]
        },

        1: {
            0: 0.5,
            1: 0.5
        },

        2: {
            0: PROBS["mutation"],
            1: 1 - PROBS["mutation"]
        }        
    }
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """

    # Identify the distribution requested
    request = dict()
    names = list(people.keys())
    for name in names:
        if name in one_gene:
            gene_no = 1
        elif name in two_genes:
            gene_no = 2
        else:
            gene_no = 0

        if name in have_trait:
            trait_prop = True
        else:
            trait_prop = False

        request[name] = {"gene":gene_no, "trait":trait_prop}

    total_prob = 1.0
    for name in names:
        current_prob = 1.0
        current_prob *= child_parent_prob(name, people, request)
        current_prob *= PROBS["trait"][request[name]["gene"]][request[name]["trait"]]

        total_prob *= current_prob


    #if request["Harry"]["trait"] == False and request["Harry"]["gene"] == 1 and request["Lily"]["trait"] == False and request["Lily"]["gene"] == 0 and request["James"]["trait"] == True and request["James"]["gene"] == 2: 
    #    print("Request of: " + str(request))
    #    print("Total probability is: " + str(total_prob))

    return total_prob

def child_parent_prob(child, people, request):
    """
    Function that returns the probability of the person having the 
    requested number of genes given their parents.
    person - the name of the person in people who the probability is being identified for
    people - the list of people with parents and traits
    request - the requested status of genes/traits
    """

    if people[child]["mother"]:
        mum = people[child]["mother"]
        dad = people[child]["father"]
        child_gene_no = request[child]["gene"]
        #mum_prob = child_parent_prob(mum, people, request)
        #dad_prob = child_parent_prob(dad, people, request)
#
        #if child_gene_no == 0:
        #    prob = mum_prob * PARENT_PROB["parent"][request[mum]["gene"]][0] * dad_prob * PARENT_PROB["parent"][request[dad]["gene"]][0]
        #elif child_gene_no == 1:
        #    prob = mum_prob * PARENT_PROB["parent"][request[mum]["gene"]][1] + dad_prob * PARENT_PROB["parent"][request[dad]["gene"]][1]
        #elif child_gene_no == 2:
        #    prob = mum_prob * PARENT_PROB["parent"][request[mum]["gene"]][1] * dad_prob * PARENT_PROB["parent"][request[dad]["gene"]][1]
        #else:
        #    raise NotImplementedError
        
        if child_gene_no == 0:
            prob = PARENT_PROB["parent"][request[mum]["gene"]][0] * PARENT_PROB["parent"][request[dad]["gene"]][0]
        elif child_gene_no == 1:
            prob = PARENT_PROB["parent"][request[mum]["gene"]][1]*(1-PARENT_PROB["parent"][request[dad]["gene"]][1]) + (1-PARENT_PROB["parent"][request[mum]["gene"]][1])*PARENT_PROB["parent"][request[dad]["gene"]][1]
        elif child_gene_no == 2:
            prob = PARENT_PROB["parent"][request[mum]["gene"]][1] * PARENT_PROB["parent"][request[dad]["gene"]][1]      

    else:

        prob = PROBS["gene"][request[child]["gene"]]

   return prob

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """

    # Identify the distribution requested
    request = dict()
    names = list(probabilities.keys())
    for name in names:
        if name in one_gene:
            gene_no = 1
        elif name in two_genes:
            gene_no = 2
        else:
            gene_no = 0

        if name in have_trait:
            trait_prop = True
        else:
            trait_prop = False

        request[name] = {"gene":gene_no, "trait":trait_prop}

    for name in names:
        probabilities[name]["gene"][request[name]["gene"]] += p
        probabilities[name]["trait"][request[name]["trait"]] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        gene_sum = probabilities[person]["gene"][0] + probabilities[person]["gene"][1] + probabilities[person]["gene"][2]
        trait_sum = probabilities[person]["trait"][False] + probabilities[person]["trait"][True]
        for gene in [0, 1, 2]:
            probabilities[person]["gene"][gene] = probabilities[person]["gene"][gene]/gene_sum

        for trait in [True, False]:
            probabilities[person]["trait"][trait] = probabilities[person]["trait"][trait]/trait_sum

if __name__ == "__main__":
    main()
