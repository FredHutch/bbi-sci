function Footer(props) {
  return React.createElement(
    "div",
    { className: "wrapper" },
    React.createElement(
      "div",
      { className: "footer-inner" },
      React.createElement(
        "div",
        { className: "footer-left" },
        React.createElement(
          "div",
          { className: "footer-logo" },
          React.createElement("img", { src: "https://brotmanbaty.org/images/BBI_Logo_REV.svg", alt: "Footer Logo" })
        )
      ),
      React.createElement(
        "div",
        { className: "footer-right" },
        React.createElement(
          "div",
          { className: "fr-wrap" },
          React.createElement(
            "p",
            null,
            "357657 | Seattle, WA 98195-8047",
            React.createElement("br", null),
            "Health Sciences Building H-564",
            React.createElement("br", null),
            "info@brotmanbaty.org",
            React.createElement("br", null),
            "206-543-9660"
          )
        )
      )
    )
  );
}

function SubFooter(props) {
  return React.createElement(
    "div",
    { className: "wrapper" },
    React.createElement(
      "span",
      { className: "subfooter-logo sl-left" },
      React.createElement("img", { src: "https://images.ctfassets.net/gad02mmq3858/7vHe9SGp6yJ3W7ZbnE1ZFH/41bf0222e2c3f0ba68f0b389f6c89243/logo_UWMedicine_Color.svg", alt: "Partner Logo" })
    ),
    React.createElement(
      "span",
      { className: "subfooter-logo sl-mid" },
      React.createElement("img", { src: "https://images.ctfassets.net/gad02mmq3858/Evsj9Rezh4aNHO4qiRtNK/a58c48d03765a90ac3bc63971e87e2a2/Group_335.svg", alt: "Partner Logo" })
    ),
    React.createElement(
      "span",
      { className: "subfooter-logo sl-right" },
      React.createElement("img", { src: "https://images.ctfassets.net/gad02mmq3858/7CSOxuV9kX7nqetiV2uXju/3ea142dec33c3a8e2649eb1c3783e26f/Group_334.svg", alt: "Partner Logo" })
    )
  );
}

ReactDOM.render(React.createElement(Footer, null), document.getElementById('footer'));

ReactDOM.render(React.createElement(SubFooter, null), document.getElementById('subfooter'));