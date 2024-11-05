import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from transformers import BertForSequenceClassification, BertTokenizerFast


def create_model():
    ckpt_path = "./checkpoints/"
    # Paths to tokenizer files
    tokenizer_files = {
        "vocab_file": ckpt_path + "vocab.txt",
        "tokenizer_file": ckpt_path + "tokenizer.json",
        "tokenizer_config_file": ckpt_path + "tokenizer_config.json",
        "special_tokens_map_file": ckpt_path + "special_tokens_map.json",
    }

    tokenizer = BertTokenizerFast(
        vocab_file=tokenizer_files["vocab_file"],
        tokenizer_file=tokenizer_files["tokenizer_file"],
        tokenizer_config=tokenizer_files["tokenizer_config_file"],
        special_tokens_map_file=tokenizer_files["special_tokens_map_file"],
    )

    model = BertForSequenceClassification.from_pretrained(ckpt_path, num_labels=10)
    model.eval()
    return model, tokenizer


def predict_and_save(model, tokenizer, smiles: str, base_save_path: str) -> np.ndarray:
    """Run a preduction and save the resulting graph and smiles image at
    - base/path/1234-graph.png
    - base/path/1234-smiles.png

    for the base path 'base/path/1234'

    """
    tokens = tokenizer(smiles, return_tensors="pt")
    predictions = model(**tokens)

    # Convert the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    mol_path = base_save_path + "-mol.png"
    # Check if the molecule was parsed correctly
    if mol is None:
        raise ValueError("Failed to parse SMILES string.")
    else:
        # Display the molecule using PIL
        # im = Draw.MolToImage(mol, size=(300, 300))
        # # im.show()

        # Alternatively, save the molecule as an image file
        Draw.MolToFile(mol, mol_path, size=(640, 480), bgColor=(0, 0, 0, 0))

    to_plot = predictions.logits.detach().numpy()[0]

    pred_file = base_save_path + "-graph.png"
    # Find the index and value of the maximum element
    max_index = np.argmax(to_plot)
    max_value = to_plot[max_index]

    # Create an array of indices for the x-axis
    indices = np.arange(len(to_plot))

    # Plot the 1D array
    plt.plot(indices, to_plot, color="blue", label="Data")

    # Highlight the maximum value
    plt.plot(max_index, max_value, "ro", label="Maximum Value")

    # Add labels and title
    plt.xlabel("Index")
    plt.ylabel("Value")
    plt.title("Logits")
    plt.legend()

    # Display the plot
    # plt.show()
    plt.savefig(pred_file, transparent=True)

    return to_plot
